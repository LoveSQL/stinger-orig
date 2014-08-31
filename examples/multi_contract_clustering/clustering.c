#include <stdint.h>

#include "stinger.h"
#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "xmalloc.h"
#include "timer.h"

double
sum_all_edgeweights(
    struct stinger * S, 
    int64_t type)
{
  double sum = 0;

  /* for each edge block */
  OMP("omp parallel for reduction(+:sum)")						
  MTA("mta assert parallel")						
  for(uint64_t eb_index = 0; eb_index < S->ETA[(type)].high; eb_index++) {	
    struct stinger_eb *  cur_eb = ebpool + S->ETA[(type)].blocks[eb_index]; 
    uint64_t stop = stinger_eb_high(cur_eb);
    for(uint64_t e = 0; e < stop; e++) { 
      if(!stinger_eb_is_blank(cur_eb, e)) {                   
	sum += cur_eb->edges[e].weight;
      }								
    } 								
  }									
}

void
sequence(uint64_t * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = i;
  }
}

void
zero(double * arr, uint64_t count) {
  OMP("omp parallel for")
  for(uint64_t i = 0; i < count; i++) {
    arr[i] = 0;
  }
}

void
score_mean_variance_first_match(
    struct stinger * S, 
    uint64_t nv, 
    double volume, 
    double * scores, 
    uint64_t * matches, 
    double * sum, 
    double * sum_squares) 
{
  double local_sum = 0;
  double local_sum_squares = 0;

  /* precompute as much of the scores as possible */
  double half_volume = volume / 2;
  double two_div_vol_sqd = 2 / ((volume) * (volume));
  struct stinger_vb * LVA = S->LVA;

  /* for each vertex */
  OMP("omp parallel for reduction(+:local_sum) reduction(+:local_sum_squares)")
  for(uint64_t u = 0; u < nv; u++) {
    double best_found = 0;
    uint64_t best_match = matches[u];

    /* precompute more of the score (vtx specific part) */
    double wt_u_x2_div_vol_sqd = LVA[u].weight * two_div_vol_sqd;

    /* for all edge blocks of that vertex */
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, u) {
      double wt_v = LVA[STINGER_EDGE_DEST].weight;
      double score = (((double)STINGER_EDGE_WEIGHT) / half_volume) - (wt_v * wt_u_x2_div_vol_sqd);
      /* check for the best score for this vertex*/
      if(score > best_found) {
	best_found = score;
	best_match = STINGER_EDGE_DEST;
      }
      /* sum the score and its square */
      local_sum += score;
      local_sum_squares += (score * score);
    } STINGER_FORALL_EDGES_OF_VTX_END();

    /* writeback the best score */
    scores[u] = best_found;
    matches[u] = best_match;
  }

  *sum = local_sum;
  *sum_squares = local_sum_squares;
}

void
filter_scores(
    uint64_t nv, 
    double * scores, 
    uint64_t * matches, 
    double sum, 
    double sum_squares) 
{
  /* compute the threshold */
  double mean = (sum) / nv;
  double mean_of_sq = (sum_squares) / nv;
  double variance = mean_of_sq - (mean * mean);
  double st_dev = sqrt(variance);
  double threshold = mean - 1.5 * st_dev;

  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    double score = scores[v];
    /* if score is below the threshold, don't match */
    if(score != 0 && score < threshold) {
      matches[v] = v;
    } else {
      /* if my match is trying to match on a lower-scoring edge, remove its match */
      uint64_t match = matches[v];
      if(scores[match] <= score) {
	matches[match] = match;
      }
    }
  }
}

void
tree_climb(
    uint64_t nv, 
    uint64_t * matches) 
{
  /* for all vertices */
  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    uint64_t older_match, old_match, match = v;
    old_match = v;
    /* climb the tree of matchings until we reach a root or cycle */
    do {
      older_match = old_match;
      old_match = match;
      match = matches[match];
      /* found a cycle - pick the lesser ID */
      if(match == older_match) {
	match = match > old_match ? old_match : match;
	break;
      }
    } while(old_match != match);
    matches[v] = match;
  }
}

int
multi_contract_root(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_vb * LVA = S->LVA;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    /* if it's a root */
    if(match_u == u) {
      /* for all edge blocks of the vertex */
      struct stinger_eb * currentBlock = ebpool + LVA[u].edges;
      while(currentBlock != ebpool) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  /* if the edge exists */
	  if(v >= 0) {
	    uint64_t match_v = matches[v];
	    /* and the edge is to be contracted, remove and increment vtx weight */
	    if(match_u == match_v) {
	      work_remaining = 1;
	      stinger_int64_fetch_add(&(LVA[match_u].weight), edge->weight);
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	      stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	    /* otherwise remove, remap, reinsert */
	    } else if(match_v != v) {
	      work_remaining = 1;
	      edge->neighbor = ~v;
	      currentBlock->numEdges--;
	      stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	      stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	      stinger_incr_edge(S, 0, u, match_v, edge->weight, timestamp);
	    }
	  }
	}
	currentBlock = currentBlock->next + ebpool;
      }
    }
  }

  return work_remaining;
}

int
multi_contract_tree(
    struct stinger * S, 
    uint64_t nv,
    uint64_t * matches,
    int64_t timestamp)
{
  int work_remaining = 0;

  struct stinger_vb * LVA = S->LVA;

  /* for all vertices */
  OMP("omp parallel for reduction(+:work_remaining)")
  for(uint64_t u = 0; u < nv; u++) {
    uint64_t match_u = matches[u];
    struct stinger_eb * currentBlock = ebpool + LVA[u].edges;
    if(match_u != u) {
      while(currentBlock != ebpool) {
	struct stinger_edge * edge = currentBlock->edges;
	struct stinger_edge * last_edge = edge + currentBlock->high;
	/* for each edge in a block */
	for(; edge < last_edge; ++edge) {
	  int64_t v = edge->neighbor;
	  if(v >= 0) {
	    work_remaining = 1;
	    uint64_t match_v = matches[v];
	    if(match_u == match_v) {
	      stinger_int64_fetch_add(&(LVA[match_u].weight), edge->weight);
	    } else {
	      stinger_incr_edge(S, 0, match_u, match_v, edge->weight, timestamp);
	    }
	    edge->neighbor = ~v;
	    currentBlock->numEdges--;
	    stinger_int64_fetch_add(&(LVA[u].outDegree), -1); 
	    stinger_int64_fetch_add(&(LVA[v].inDegree), -1);
	  }
	}
	currentBlock = currentBlock->next + ebpool;
      }
      stinger_int64_fetch_add(&S->LVA[match_u].weight, S->LVA[u].weight);
      S->LVA[u].weight = 0;
    }
  }

  return work_remaining;
}

double
modularity_score (
    struct stinger * S, 
    uint64_t * cluster, 
    uint64_t nv, 
    uint64_t ne)
{
  int64_t *Lss = (int64_t *)calloc(nv, sizeof(int64_t)), *Lsplus = (int64_t *)calloc(nv, sizeof(int64_t));
  double mod = 0.0, m_inv = 1.0/(double)ne;

  STINGER_FORALL_EDGES_BEGIN(S, 0) {
    if(STINGER_EDGE_SOURCE < STINGER_EDGE_DEST) {
      if(cluster[STINGER_EDGE_SOURCE] == cluster[STINGER_EDGE_DEST]) {
	Lss[cluster[STINGER_EDGE_SOURCE]]++;
      } else {
	Lsplus[cluster[STINGER_EDGE_SOURCE]]++;
      }
    }
  } STINGER_FORALL_EDGES_END();

  uint64_t total = 0;
  for(uint64_t i = 0; i < nv; i++) {
    if(cluster[i] >= nv)
      fprintf(stderr,"\nERROR: CLUSTER of %ld is %ld",i,cluster[i]);
    if(cluster[i] == i)
      total++;
    if(Lss[i] != 0.0 || Lsplus[i] != 0.0) {
      mod += ((double)Lss[i]) - (m_inv * ((double)Lsplus[i])) * (0.25 * ((double)Lsplus[i]));
    }
  }

  free(Lss); free(Lsplus);

  return mod * m_inv;
}

void
static_multi_contract_clustering (
    uint64_t ** matches,
    uint64_t nv,
    struct stinger * S)
{

  double volume = 0;

  volume = sum_all_edgeweights(S, 0);

  int work_remaining = 1;

  uint64_t iteration = 0;
  double sum, sum_squares;
  double * scores = xmalloc(nv * sizeof(double));
  *matches = xmalloc(nv * sizeof(uint64_t));

  OMP("omp parallel for")
  for(uint64_t v = 0; v < nv; v++) {
    stinger_set_vweight(S,v,0);
  }

  sequence(*matches, nv);

  int64_t max_iterations = 1000;
  if(getenv("MAX_ITERATIONS") && atoi(getenv("MAX_ITERATIONS")) > 0)
     max_iterations = atoi(getenv("MAX_ITERATIONS")) - 1;

  PRINT_STAT_INT64("max_iterations", max_iterations);

  tic();
  while(work_remaining && iteration < max_iterations) {
    iteration++;
    //PRINT_STAT_INT64("starting_iteration", iteration); fflush(stdout);

    zero(scores, nv);
    work_remaining = 0;
    sum = 0;
    sum_squares = 0;

    score_mean_variance_first_match(S, nv, volume, scores, *matches, &sum, &sum_squares);

    filter_scores(nv, scores, *matches, sum, sum_squares);

    tree_climb(nv, *matches); 

    work_remaining += multi_contract_root(S, nv, *matches, 0);

    work_remaining += multi_contract_tree(S, nv, *matches, 0);
  }
  PRINT_STAT_DOUBLE("cluster_time", toc());

  uint64_t count = 0;
  uint64_t edges = 0;
  for(uint64_t v = 0; v < nv; v++) {
    if((*matches)[v] == v)
      count++;
    if(stinger_outdegree(S, v))
      edges += stinger_outdegree(S, v);
  }

  PRINT_STAT_INT64("clusters", count);
  PRINT_STAT_INT64("clusters_edges", edges);

  PRINT_STAT_INT64("clustering_iterations", iteration);
}


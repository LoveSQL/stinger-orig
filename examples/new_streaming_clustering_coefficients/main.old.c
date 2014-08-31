/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger.h"
#include "timer.h"
#include "xmalloc.h"

#define ACTI(k) (action[2*(k)])
#define ACTJ(k) (action[2*(k)+1])

static int64_t nv, ne, naction;
static int64_t * restrict off;
static int64_t * restrict from;
static int64_t * restrict ind;
static int64_t * restrict weight;
static int64_t * restrict action;

/* handles for I/O memory */
static int64_t * restrict graphmem;
static int64_t * restrict actionmem;

static char * initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
static char * action_stream_name = ACTION_STREAM_NAME_DEFAULT;

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;

static struct stinger * S;

void count_all_triangles (const size_t nv, const int64_t * off, const int64_t * ind, int64_t * ntri);

static int
i64cmp (const void *ap, const void *bp)
{
  return (int) ((*(int64_t *) ap) - (*(int64_t *) bp));
}

int
main (const int argc, char *argv[])
{
  parse_args (argc, argv, &initial_graph_name, &action_stream_name, &batch_size, &nbatch);
  STATS_INIT();

  load_graph_and_action_stream (initial_graph_name, &nv, &ne, (int64_t**)&off,
	      (int64_t**)&ind, (int64_t**)&weight, (int64_t**)&graphmem,
	      action_stream_name, &naction, (int64_t**)&action, (int64_t**)&actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, -2);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);

  int64_t * ntri = xcalloc (nv, sizeof (*ntri));
  int64_t * local_cc = xcalloc (nv, sizeof (*local_cc));
  int64_t * affected = xcalloc (nv, sizeof (*affected));

  /* Calculate the initial clustering coefficients */
  tic();
  count_all_triangles (nv, off, ind, ntri);
  PRINT_STAT_DOUBLE ("time_triangles", toc());
  fflush(stdout);

  tic();
  int64_t global_degsum = 0;
  int64_t global_ntri = 0;
  OMP("omp parallel for")
  MTA("mta assert nodep")
  for (int64_t i = 0; i < nv; i++) {
    int64_t deg = off[i+1] - off[i];
    int64_t d = deg * (deg-1);
    local_cc[i] = (d ? ntri[i] / (double) d : 0.0);
    global_degsum += d;
    global_ntri += ntri[i];
  }
  double global_cc = (global_degsum ? global_ntri / (double) global_degsum : 0.0);
  PRINT_STAT_DOUBLE ("initial_cc", global_cc);
  PRINT_STAT_DOUBLE ("time_initial_cc", toc());
  PRINT_STAT_HEX64("error_code", stinger_consistency_check (S, nv));
  fflush(stdout);

  double *update_time_trace = xmalloc (nbatch * sizeof(*update_time_trace));

  /* Updates */
  int64_t ntrace = 0;
  if (batch_size == 1)
  {
    for (int64_t actk = 0; actk < nbatch; actk++, ntrace++) {
      tic();

      int64_t i = ACTI(actk);
      int64_t j = ACTJ(actk);

      int64_t changed;
      int64_t dir;
      int64_t incr = 0;
      int64_t num_affected = 0;
      int64_t delta_global_ntri;

      if (i < 0) {  /* This is a delete */
	i = -i-1;
	j = -j-1;
	dir = -1;

	assert (i >= 0);
	assert (i < nv);
	assert (j >= 0);
	assert (j < nv);

	changed = stinger_remove_edge_pair (S, 0, i, j);

      } else {  /* This is an insert */
	dir = 1;

	assert (i >= 0);
	assert (i < nv);
	assert (j >= 0);
	assert (j < nv);

	changed = stinger_insert_edge_pair (S, 0, i, j, 1, ntrace+1);

      }

      if (changed && i != j) {
	affected[0] = i;
	affected[1] = j;
	num_affected = 2;

	int64_t i_deg, j_deg;
	i_deg = stinger_outdegree(S, i);
	j_deg = stinger_outdegree(S, j);

	if (j_deg && i_deg >= j_deg) {
	  STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, i) {
	    int64_t neighbor_i = STINGER_EDGE_DEST;
	    if (neighbor_i != i) {
	      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, j) {
		int64_t neighbor_j = STINGER_EDGE_DEST;
		if (j != neighbor_j && neighbor_i == neighbor_j) {
		  incr += dir;
		  stinger_int64_fetch_add (&ntri[neighbor_i], 2*dir);
		  affected[stinger_int64_fetch_add (&num_affected, 1)] = neighbor_i;
#if !defined(__MTA__)
		  break;
#endif
		}
	      } STINGER_FORALL_EDGES_OF_VTX_END();
	    }
	  } STINGER_FORALL_EDGES_OF_VTX_END();
	} else if (i_deg) {
	  STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, j) {
	    int64_t neighbor_j = STINGER_EDGE_DEST;
	    if (neighbor_j != j) {
	      STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, i) {
		int64_t neighbor_i = STINGER_EDGE_DEST;
		if (i != neighbor_i && neighbor_j == neighbor_i) {
		  incr += dir;
		  stinger_int64_fetch_add (&ntri[neighbor_j], 2*dir);
		  affected[stinger_int64_fetch_add (&num_affected, 1)] = neighbor_j;
#if !defined(__MTA__)
		  break;
#endif
		}
	      } STINGER_FORALL_EDGES_OF_VTX_END();
	    }
	  } STINGER_FORALL_EDGES_OF_VTX_END();
	}

	ntri[i] += 2*incr;
	ntri[j] += 2*incr;
	delta_global_ntri = 6*incr;

	/* Update the local clustering coefficient for affected vertices */
	for (int64_t k = 0; k < num_affected; k++) {
	  int64_t v = affected[k];
	  double d = (double) stinger_outdegree (S, v);
	  double new_local_cc = (d > 1 ? ntri[v] / (d * (d-1)) : 0.0);
	  local_cc[v] = new_local_cc;  // XXX is this safe?
	}

	/* Update the global clustering coefficient */
	i_deg = stinger_outdegree (S, i);
	j_deg = stinger_outdegree (S, j);

	int64_t gntri = stinger_int64_fetch_add (&global_ntri, delta_global_ntri) + delta_global_ntri;
	incr = 2 * (i_deg + j_deg - 2);
	if (dir < 0)  // removed an edge
	  incr = -2 * (i_deg + j_deg);
	else  // inserted an edge
	  incr = 2 * (i_deg + j_deg - 2);
	int64_t gdegsum = stinger_int64_fetch_add (&global_degsum, incr) + incr;
	global_cc = (gdegsum ? gntri / (double) gdegsum : 0.0);
      }

      update_time_trace[ntrace] = toc();
    }

  } else { /* batching */

    int64_t * act = xmalloc (2 * 2 * batch_size * sizeof(*act));
    int64_t * insoff = xmalloc ((2 * 2 * batch_size + 1) * sizeof(*insoff));
    int64_t * deloff = xmalloc ((2 * 2 * batch_size + 1) * sizeof(*deloff));
    int64_t * result = xmalloc (batch_size * sizeof(*result));

    for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
    {
      tic();

      int64_t delta_global_ntri = 0;
      int64_t num_affected = 0;

      const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
      int64_t N = stinger_sort_actions (endact - actno, &action[2*actno], insoff, deloff, act);

      OMP("omp parallel for")
      MTA("mta assert parallel")
      for (int64_t k = 0; k < nv; k++)
      {
	affected[k] = 0;
      }

      OMP("omp parallel for")
      MTA("mta assert parallel")
      for (int64_t k = 0; k < batch_size; k++)
      {
	result[k] = 0;
      }

      /* This section processes all of the pending deletes to determine which triangles they will break.
	 We are still working on a way to efficiently detect multiple incident deletes so that they are not
	 double counted. */

      OMP("omp parallel for")
      MTA("mta assert nodep")
      for (int64_t k = 0; k < N; k++)
      {
	int64_t i = act[2*deloff[k]];

	int64_t i_deg, j_deg;
	i_deg = stinger_outdegree(S, i);
	
	int64_t * i_neighbors = xmalloc (i_deg * sizeof(*i_neighbors));
	size_t i_degree;
	stinger_gather_typed_successors (S, 0, i, &i_degree, i_neighbors, i_deg);
	assert (i_degree == i_deg);
	qsort (i_neighbors, i_degree, sizeof(i_neighbors[0]), i64cmp);

	for (int64_t w = deloff[k]; w < insoff[k]; w++)
	{
	  int64_t incr = 0;
	  int64_t dir = -1;

	  if (w > deloff[k] && act[2*w+1] == act[2*(w-1)+1]) continue;  // don't process multiple copies of the same edge action

	  int64_t j = act[2*w+1];

	  if (i > j) continue;

	  if (find_in_sorted (j, i_degree, i_neighbors) == -1) continue;
	  //if (stinger_has_typed_successor (S, 0, i, j) == 0) continue;  // don't process if it's not in the graph

	  stinger_int64_fetch_add (&affected[i], 1);
	  stinger_int64_fetch_add (&affected[j], 1);

	  j_deg = stinger_outdegree(S, j);
	  size_t j_degree;

	  if (!j_deg)
	    continue;

	  int64_t * j_neighbors = xmalloc (j_deg * sizeof(*j_neighbors));
	  stinger_gather_typed_successors (S, 0, j, &j_degree, j_neighbors, j_deg);
	  assert (j_degree == j_deg);

	  for (int64_t k2 = 0; k2 < j_degree; k2++)
	  {
	    int64_t u = j_neighbors[k2];
	    int64_t where = find_in_sorted (u, i_degree, i_neighbors);
	    if (where >= 0) {
	      incr += dir;
	      stinger_int64_fetch_add (&ntri[u], 2 * dir);
	      stinger_int64_fetch_add (&affected[u], 1);
	    }
	  }

	  free (j_neighbors);
	}

	free (i_neighbors);
      }

      /* For each unique vertex represented in the batch, we remove its degree contribution
	 from the global degree sum, or the denominator of the global clustering coefficient. */

      OMP("omp parallel for")
      MTA("mta assert nodep")
      for (int64_t k = 0; k < N; k++)
      {
	if (deloff[k] < insoff[k]) {
	  int64_t i = act[2*deloff[k]];
	  int64_t deg = stinger_outdegree (S, i);
	  stinger_int64_fetch_add (&global_degsum, -deg*(deg-1));
	}
      }

      int64_t *actions = &action[2*actno];
      int64_t numActions = endact - actno;

      MTA("mta assert parallel")
      MTA("mta block dynamic schedule")
      OMP("omp parallel for")
      for(uint64_t k = 0; k < endact - actno; k++) {
	const int64_t i = actions[2 * k];
	const int64_t j = actions[2 * k + 1];

	if (i != j && i < 0) {
	  int64_t r1 = stinger_remove_edge(S, 0, ~i, ~j);
	  int64_t r2 = stinger_remove_edge(S, 0, ~j, ~i);
	}

	if (i != j && i >= 0) {
	  int64_t r1 = stinger_insert_edge (S, 0, i, j, 1, actno+1);
	  int64_t r2 = stinger_insert_edge (S, 0, j, i, 1, actno+1);

	  if (r1==1 && r2==1) result[k] = 3;
	  else if (r1==1) result[k] = 1;
	  else if (r2==1) result[k] = 2;
	}
      }

      /* This appears to be broken */ 
      //stinger_remove_and_insert_batch (S, 0, actno+1, N, insoff, deloff, act);

      /* Now that insertions and deletions are processed, we update the global degree sum with
	 the new degree counts. */

      OMP("omp parallel for")
      MTA("mta assert nodep")
      for (int64_t k = 0; k < N; k++)
	if (deloff[k] < insoff[k]) {
	  int64_t i = act[2*deloff[k]];
	  int64_t deg = stinger_outdegree (S, i);
	  stinger_int64_fetch_add (&global_degsum, +deg*(deg-1));
	}

      /* In this step, we will process the insertions to determine how many triangles were
	 created.  Multiple copies of the same insertion are ignored.  An ordering is
	 determined to ensure that multiple insertions incident on the same vertex are
	 considered properly.

	 Note that some insertions may have already been in the graph, and we use the first
	 timestep to determine if that's the case. */

      MTA("mta assert parallel")
      MTA("mta block dynamic schedule")
      for(uint64_t k = 0; k < endact - actno; k++) {
	const int64_t i = actions[2 * k];
	const int64_t j = actions[2 * k + 1];

	if (i != j && i >= 0) {
	  int64_t i_deg, j_deg;
	  i_deg = stinger_outdegree(S, i);

	  int64_t * i_neighbors = xmalloc (i_deg * sizeof(*i_neighbors));
	  size_t i_degree;
	  stinger_gather_typed_successors (S, 0, i, &i_degree, i_neighbors, i_deg);
	  assert (i_degree == i_deg);
	  qsort (i_neighbors, i_degree, sizeof(i_neighbors[0]), i64cmp);

	  int64_t incr = 0;

	  if (!result[k]) continue;

	  stinger_int64_fetch_add (&affected[i], 1);
	  stinger_int64_fetch_add (&affected[j], 1);
	  j_deg = stinger_outdegree(S, j);
	  size_t j_degree;

	  if (!j_deg)
	    continue;

	  int64_t * j_neighbors = xmalloc (j_deg * sizeof(*j_neighbors));
	  stinger_gather_typed_successors (S, 0, j, &j_degree, j_neighbors, j_deg);
	  assert (j_degree == j_deg);

	  for (int64_t k2 = 0; k2 < j_degree; k2++)
	  {
	    int64_t u = j_neighbors[k2];
	    int64_t where = find_in_sorted (u, i_degree, i_neighbors);
	    if (where >= 0) {
	      int64_t neighbor_i_ts = stinger_edge_timestamp_first (S, i, u, 0);
	      int64_t neighbor_j_ts = stinger_edge_timestamp_first (S, j, u, 0);
	      if (neighbor_i_ts != actno+1 && neighbor_j_ts != actno+1) {
		incr += 1;
		stinger_int64_fetch_add (&ntri[u], 2);
		stinger_int64_fetch_add (&affected[u], 1);
	      } else if (i < u) {
		if (j < u || neighbor_i_ts != actno+1) {
		  incr += 1;
		  stinger_int64_fetch_add (&ntri[u], 2);
		  stinger_int64_fetch_add (&affected[u], 1);
		} 
	      }
	    }
	  }

	  free (j_neighbors);

	  stinger_int64_fetch_add (&ntri[i], 2*incr);
	  stinger_int64_fetch_add (&ntri[j], 2*incr);
	  delta_global_ntri += (6*incr);

	  free (i_neighbors);
	}
      }

      /* Update the local clustering coefficient for affected vertices */
      OMP("omp parallel for")
      MTA("mta assert parallel")
      for (int64_t k = 0; k < nv; k++) {
	if (affected[k]) {
	  double d = (double) stinger_outdegree (S, k);
	  double new_local_cc = (d > 1 ? ntri[k] / (d * (d-1)) : 0.0);
	  local_cc[k] = new_local_cc;
	}
      }

      global_ntri += delta_global_ntri;
      update_time_trace[ntrace] = toc();
      ++ntrace;
    } /* End of batch */

    free (deloff);
    free (insoff);
    free (act);

  }

  /* Print the times */
  double time_updates = 0;
  for (int64_t k = 0; k < nbatch; k++) {
    time_updates += update_time_trace[k];
  }
  global_cc = (global_degsum ? global_ntri / (double) global_degsum : 0.0);
  PRINT_STAT_DOUBLE ("final_cc", global_cc);
  PRINT_STAT_DOUBLE ("time_updates", time_updates);
  PRINT_STAT_DOUBLE ("updates_per_sec", (nbatch * batch_size) / time_updates); 

  PRINT_STAT_INT64 ("global_ntri", global_ntri);
  PRINT_STAT_INT64 ("global_degsum", global_degsum);

  PRINT_STAT_HEX64("error_code", stinger_consistency_check (S, nv));

  free (update_time_trace);
  free (affected);
  free (local_cc);
  free (ntri);

  stinger_free_all (S);
  free (graphmem);
  free (actionmem);
  STATS_END();
}



MTA("mta expect parallel context") MTA("mta inline")
size_t
count_intersections (const int64_t ai,
		     const size_t alen,
		     const int64_t * a,
		     const int64_t bi,
		     const size_t blen,
		     const int64_t * b)
{
  size_t ka = 0, kb = 0;
  size_t out = 0;

  if (!alen || !blen || a[alen-1] < b[0] || b[blen-1] < a[0])
    return 0;

  while (1) {
    if (ka >= alen || kb >= blen) break;

    int64_t va = a[ka];
    int64_t vb = b[kb];

    /* Skip self-edges. */
    if (UNLIKELY(va == ai)) {
      ++ka;
      if (ka >= alen) break;
      va = a[ka];
    }
    if (UNLIKELY(vb == bi)) {
      ++kb;
      if (kb >= blen) break;
      vb = b[kb];
    }

    if (va == vb) { ++ka; ++kb; ++out; }
    else if (va < vb) {
      /* Advance ka */
      ++ka;
      while (ka < alen && a[ka] < vb) ++ka;
    } else {
      /* Advance kb */
      ++kb;
      while (kb < blen && va > b[kb]) ++kb;
    }
  }

  return out;
}


MTA ("mta expect parallel context") MTA ("mta inline") int64_t
count_triangles (const size_t nv, const int64_t * restrict off,
                 const int64_t * restrict ind, int64_t i)
/* Assume sorted index lists. */
{
  int64_t out = 0;
  size_t k;
  const int64_t *restrict iind = &ind[off[i]];
  const size_t ideg = off[i + 1] - off[i];
  const size_t k_end = off[i + 1];

  MTA ("mta assert nodep") MTA ("mta assert noalias")
    for (size_t k = 0; k < ideg; ++k) {
      const int64_t j = iind[k];
      if (j != i) {
        const int64_t *restrict jind = &ind[off[j]];
        const size_t jdeg = off[j + 1] - off[j];
        const size_t shared =
          count_intersections (i, ideg, iind, j, jdeg, jind);
        out += (int64_t) shared;
      }
    }
  return out;
}


MTA ("mta expect serial context")
void
count_all_triangles (const size_t nv, const int64_t * restrict off,
                     const int64_t * restrict ind,
                     int64_t * restrict ntri)
{
  {
    const size_t N = nv;
    OMP ("omp for schedule(dynamic,128)")
    MTA ("mta dynamic schedule") MTASTREAMS ()
    for (size_t i = 0; i < N; ++i)
      ntri[i] = count_triangles (nv, off, ind, i);
  }
}

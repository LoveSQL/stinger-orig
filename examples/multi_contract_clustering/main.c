/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger.h"
#include "timer.h"
#include "xmalloc.h"
#include "clustering.h"
#include <omp.h>

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
static struct stinger * S_cluster;

static double * update_time_trace;

int64_t parallel_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                              int64_t * component_map);

void make_clique(struct stinger * S, uint64_t start, uint64_t size) {
  for(uint64_t i = start; i < start+size; i++) {
    for(uint64_t j = i+1; j < start+size; j++) {
      stinger_insert_edge(S, 0, i, j, 1, 0);
      stinger_insert_edge(S, 0, j, i, 1, 0);
    }
  }
}

void make_clique_ring(struct stinger * S, uint64_t start, uint64_t clique_size, uint64_t num_cliques) {
  for(uint64_t i = 0; i < num_cliques; i++) {
    make_clique(S, start+i*clique_size, clique_size);
  }
  for(uint64_t i = 0; i < num_cliques - 1; i++) {
    stinger_insert_edge(S, 0, start + i*clique_size, start + (i+1)*clique_size, 1, 0);
    stinger_insert_edge(S, 0, start + (i+1)*clique_size, start + i*clique_size, 1, 0);
  }
  stinger_insert_edge(S, 0, start, start + (num_cliques-1)*clique_size, 1, 0);
  stinger_insert_edge(S, 0, start + (num_cliques-1)*clique_size, start, 1, 0);
}

#define LINE_SIZE 1048576
void load_text_adjacency(struct stinger * S, const char * file_name, uint64_t * nv, uint64_t * ne) {
  FILE * fp = fopen(file_name, "r");
  int64_t * off;
  int64_t * ind;
  int64_t * weight;
  char line[LINE_SIZE];

  if(fp) {
    /* eat comments */
    line[0] = '%';
    while(line[0] == '%') {
      if(!fgets(line, LINE_SIZE, fp))
	return;
    }
    sscanf(line, "%ld %ld", nv, ne);

    for(uint64_t u = 0; fgets(line, LINE_SIZE, fp); u++) {
      /* eat comments */
      if(line[0] == '%') {
	u--;
	continue;
      }

      uint64_t v = 0;
      char * ptr = line;
      int read = 0;
      while(sscanf(ptr, "%ld%n", &v, &read) > 0) {
	ptr += read;
	stinger_insert_edge(S, 0, u, v-1, 1, 0);
      }
    }

    fclose(fp);
  }
}

int
main (const int argc, char *argv[])
{
  STATS_INIT();
#if 0
  parse_args (argc, argv, &initial_graph_name, &action_stream_name, &batch_size, &nbatch);

  load_graph_and_action_stream (initial_graph_name, &nv, &ne, (int64_t**)&off,
	      (int64_t**)&ind, (int64_t**)&weight, (int64_t**)&graphmem,
	      action_stream_name, &naction, (int64_t**)&action, (int64_t**)&actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK();

  update_time_trace = xmalloc (nbatch * sizeof(*update_time_trace));

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  S_cluster = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  stinger_set_initial_edges (S_cluster, nv, 0, off, ind, weight, NULL, NULL, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);
  free(graphmem);
#elif 0
  tic ();
  S = stinger_new ();
  S_cluster = stinger_new ();
  uint64_t p_size = atoi(getenv("PROBLEM_SIZE"));
  make_clique(S, 0, p_size);
  make_clique(S, p_size, p_size);
  stinger_insert_edge(S, 0, p_size-1, p_size, 1, 0);
  stinger_insert_edge(S, 0, p_size, p_size-1, 1, 0);
  make_clique(S_cluster, 0, p_size);
  make_clique(S_cluster, p_size, p_size);
  stinger_insert_edge(S_cluster, 0, p_size-1, p_size, 1, 0);
  stinger_insert_edge(S_cluster, 0, p_size, p_size-1, 1, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  nv = p_size*2;
#elif 0
  tic ();
  S = stinger_new ();
  S_cluster = stinger_new ();
  uint64_t c_size = atoi(getenv("CLIQUE_SIZE"));
  uint64_t nc = atoi(getenv("NUM_CLIQUES"));
  make_clique_ring(S, 0, c_size, nc);
  make_clique_ring(S_cluster, 0, c_size, nc);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  nv = nc * c_size;
#else 
  S = stinger_new();
  load_text_adjacency(S, getenv("GRAPH_FILE"), &nv, &ne);
  S_cluster = stinger_new();
  load_text_adjacency(S_cluster, getenv("GRAPH_FILE"), &nv, &ne);
  printf("\n\t\"graph_name\": \"%s\"", getenv("GRAPH_FILE"));
#endif

  OMP("omp parallel")
  {
    OMP("omp master")
    {
      PRINT_STAT_INT64("num_threads", (int64_t)omp_get_num_threads());
    }
  }


  tic ();
  uint64_t errorCode = stinger_consistency_check (S, nv);
  double time_check = toc ();
  PRINT_STAT_HEX64("error_code", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  /* Serial connected components */
  uint64_t * component_map = xmalloc(sizeof(uint64_t) * nv);
  tic ();
  int64_t num_comp =
    parallel_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_parallel_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);

  /* Static clustering */
  uint64_t * matches;
  static_multi_contract_clustering(&matches, nv, S_cluster);
  PRINT_STAT_DOUBLE("clustering_modularity", modularity_score(S, matches, nv, stinger_total_edges(S)));

  /* Print the times */
  double time_updates = 0;
  for (int64_t k = 0; k < nbatch; k++) {
    time_updates += update_time_trace[k];
  }
  //PRINT_STAT_DOUBLE ("time_updates", time_updates);
  //PRINT_STAT_DOUBLE ("updates_per_sec", (nbatch * batch_size) / time_updates); 

  tic ();
  errorCode = stinger_consistency_check (S, nv);
  time_check = toc ();
  PRINT_STAT_HEX64("error_code", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  PRINT_STAT_INT64_ARRAY_AS_PAIR("cluster_mapping", matches, nv);

  free(update_time_trace);
  stinger_free_all (S);
  free (matches);
  free (actionmem);
  STATS_END();
}

/*
 * Perform a shiloach vishkin connected components calculation in parallel on a stinger graph
 * and return the number of components found
 */
int64_t
parallel_shiloach_vishkin_components (struct stinger * S, int64_t nv,
                                      int64_t * component_map)
{
  /* Initialize each vertex with its own component label in parallel */
  OMP ("omp parallel for")
    for (uint64_t i = 0; i < nv; i++) {
      component_map[i] = i;
    }

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

    /* For all edges in the STINGER graph of type 0 in parallel, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
    STINGER_PARALLEL_FORALL_EDGES_BEGIN (S, 0) {
      if (component_map[STINGER_EDGE_DEST] <
          component_map[STINGER_EDGE_SOURCE]) {
        component_map[STINGER_EDGE_SOURCE] = component_map[STINGER_EDGE_DEST];
        changed++;
      }
    }
    STINGER_PARALLEL_FORALL_EDGES_END ();

    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing with OpenMP parallel for */
    OMP ("omp parallel for")
      MTA ("mta assert nodep")
      for (uint64_t i = 0; i < nv; i++) {
        while (component_map[i] != component_map[component_map[i]])
          component_map[i] = component_map[component_map[i]];
      }
  }

  /* Count components */
  uint64_t components = 1;
  MTA ("mta assert nodep")
    OMP ("omp parallel for reduction(+:components)")
    for (uint64_t i = 1; i < nv; i++) {
      if (component_map[i] == i) {
        components++;
      }
    }

  return components;
}

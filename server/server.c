/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger-shared.h"
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

static double * update_time_trace;

int64_t
parallel_shiloach_vishkin_components (struct stinger * S, int64_t nv,
                                      int64_t * component_map);

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

  update_time_trace = xmalloc (nbatch * sizeof(*update_time_trace));

  /* Convert to STINGER */
  tic ();
  struct stinger_shared * shared;
  char * name;
  S = stinger_shared_new (&shared, &name);
  printf("\n\tNAME: %s", name);
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);

  free(graphmem);

  tic ();
  uint32_t errorCode = stinger_consistency_check (S, nv);
  double time_check = toc ();
  printf("\n\t\"error_code\": 0x%lx", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  /* Updates */
  int64_t ntrace = 0;

  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
  {
    tic();

    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);
    int64_t *actionStream = &action[2*actno];
    int64_t numActions = endact - actno;

    for (int64_t k = 0; k < numActions; k++) {
      const int64_t i = actionStream[2 * k];
      const int64_t j = actionStream[2 * k + 1];

      /* This is where you work on each edge action. 
       * If (i >= 0 && j >= 0) then (i, j) is an edge insertion.
       * If (i < 0 ** j < 0) then (~i, ~j) is an edge deletion.
       */

    }

    update_time_trace[ntrace] = toc();
    ntrace++;

  } /* End of batch */

  int64_t *component_map = xmalloc (nv * sizeof (int64_t));
  tic ();
  int64_t num_comp = parallel_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_parallel_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);
  free(component_map);

  getchar();



  /* Print the times */
  double time_updates = 0;
  for (int64_t k = 0; k < nbatch; k++) {
    time_updates += update_time_trace[k];
  }
  PRINT_STAT_DOUBLE ("time_updates", time_updates);
  PRINT_STAT_DOUBLE ("updates_per_sec", (nbatch * batch_size) / time_updates); 

  tic ();
  errorCode = stinger_consistency_check (S, nv);
  time_check = toc ();
  printf("\n\t\"error_code\": 0x%lx", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  free(update_time_trace);

  stinger_shared_free (S, shared, name);
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

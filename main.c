/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger.h"
#include "stinger-iterator.h"
#include "timer.h"
#include "xmalloc.h"

#define ACTI(k) (action[2*(k)])
#define ACTJ(k) (action[2*(k)+1])

static int64_t nv, ne, naction;
static int64_t *restrict off;
static int64_t *restrict from;
static int64_t *restrict ind;
static int64_t *restrict weight;
static int64_t *restrict action;

/* handles for I/O memory */
static int64_t *restrict graphmem;
static int64_t *restrict actionmem;

static char *initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
static char *action_stream_name = ACTION_STREAM_NAME_DEFAULT;

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;

static struct stinger *S;

int64_t serial_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                            int64_t * component_map);
int64_t parallel_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                              int64_t * component_map);
int64_t
serial_iterator_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                    int64_t * component_map);

int
main (const int argc, char *argv[])
{
  parse_args (argc, argv, &initial_graph_name, &action_stream_name,
              &batch_size, &nbatch);
  STATS_INIT ();

  load_graph_and_action_stream (initial_graph_name, &nv, &ne,
                                (int64_t **) & off, (int64_t **) & ind,
                                (int64_t **) & weight,
                                (int64_t **) & graphmem, action_stream_name,
                                &naction, (int64_t **) & action,
                                (int64_t **) & actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK ();

  int64_t *component_map = xmalloc (nv * sizeof (int64_t));

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush (stdout);

  /* Serial connected components */
  tic ();
  int64_t num_comp =
    serial_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_serial_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);

  tic ();
  num_comp =
    serial_iterator_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_serial_iterator_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);

  /* Parallel connected components */
  tic ();
  num_comp = parallel_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_parallel_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);

  stinger_free_all (S);
  free (graphmem);
  free (actionmem);
  STATS_END ();
}

/*
 * Perform a shiloach vishkin connected components calculation on a stinger graph
 * and return the number of components found
 */
int64_t
serial_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                    int64_t * component_map)
{
  /* Initialize each vertex with its own component label in parallel */
  for (uint64_t i = 0; i < nv; i++) {
    component_map[i] = i;
  }

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

    /* For all edges in the STINGER graph of type 0 in parallel, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
    STINGER_FORALL_EDGES_BEGIN (S, 0) {
      if (component_map[STINGER_EDGE_DEST] <
          component_map[STINGER_EDGE_SOURCE]) {
        component_map[STINGER_EDGE_SOURCE] = component_map[STINGER_EDGE_DEST];
        changed++;
      }
    }
    STINGER_FORALL_EDGES_END ();

    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing with OpenMP parallel for */
    for (uint64_t i = 0; i < nv; i++) {
      while (component_map[i] != component_map[component_map[i]])
        component_map[i] = component_map[component_map[i]];
    }
  }

  /* Count components */
  uint64_t components = 1;
  for (uint64_t i = 1; i < nv; i++) {
    if (component_map[i] == i) {
      components++;
    }
  }

  return components;
}

/*
 * Perform a shiloach vishkin connected components calculation on a stinger graph
 * and return the number of components found
 */
int64_t
serial_iterator_shiloach_vishkin_components (struct stinger *S, int64_t nv,
                                    int64_t * component_map)
{
  /* Initialize each vertex with its own component label */
  for (uint64_t i = 0; i < nv; i++) {
    component_map[i] = i;
  }

  /* Initialize an iterator */
  int64_t types = 0, count = 1;
  stinger_iterator_t * iter = stinger_iterator_new(S);
  stinger_iterator_edge_type_filter_no_copy(&types, count, iter);

  /* Iterate until no changes occur */
  while (1) {
    int changed = 0;

    /* For all edges in the STINGER graph of type 0, attempt to assign
       lesser component IDs to neighbors with greater component IDs */
    while(stinger_iterator_next(iter)) {
      if (component_map[iter->dest] < component_map[iter->source]) {
        component_map[iter->source] = component_map[iter->dest];
        changed++;
      }
    }

    /* if nothing changed */
    if (!changed)
      break;

    /* Tree climbing */
    for (uint64_t i = 0; i < nv; i++) {
      while (component_map[i] != component_map[component_map[i]])
        component_map[i] = component_map[component_map[i]];
    }
  }

  /* Count components */
  uint64_t components = 1;
  for (uint64_t i = 1; i < nv; i++) {
    if (component_map[i] == i) {
      components++;
    }
  }

  return components;
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

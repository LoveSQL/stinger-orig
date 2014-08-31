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

int64_t
parallel_shiloach_vishkin_components (struct stinger * S, int64_t nv,
                                      int64_t * component_map);

static struct stinger * S;

int
main (const int argc, char *argv[])
{

  printf("\n\tMapping: %s", argv[1]);
  /* Convert to STINGER */
  struct stinger_shared * shared;
  S = stinger_shared_map (&shared, argv[1]);
  uint64_t nv = stinger_max_active_vertex(S) + 1;


  tic ();
  uint32_t errorCode = stinger_consistency_check (S, nv);
  double time_check = toc ();
  printf("\n\t\"error_code\": 0x%lx", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  int64_t *component_map = xmalloc (nv * sizeof (int64_t));
  tic ();
  int64_t num_comp = parallel_shiloach_vishkin_components (S, nv, component_map);
  PRINT_STAT_DOUBLE ("time_parallel_components", toc ());
  PRINT_STAT_INT64 ("number_of_components", num_comp);
  free(component_map);

  tic ();
  errorCode = stinger_consistency_check (S, atoi(argv[1]));
  time_check = toc ();
  printf("\n\t\"error_code\": 0x%lx", errorCode);
  PRINT_STAT_DOUBLE ("time_check", time_check);

  stinger_shared_unmap (S, shared, argv[1]);
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

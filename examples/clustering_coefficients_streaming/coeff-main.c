/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "stinger.h"
#include "timer.h"
#include "xmalloc.h"
#include "clustering-coeff.h"

#define ACTI(k) (action[2*(k)])
#define ACTJ(k) (action[2*(k)+1])
#define ACTI2(k) (act[2*(k)])
#define ACTJ2(k) (act[2*(k)+1])

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

#define MAX_DISAGREE_PRINTED 10


extern void simple_update_tris (int dir, struct stinger *G,
                                int64_t v, int64_t w,
                                int64_t * naffected_in,
                                int64_t * affected,
                                int64_t * ntri,
                                int64_t * delta_global_ntri_in,
                                int64_t * gather_work);
extern void bulk_update_tris (struct stinger *G, int64_t v,
                              int64_t nedge, int64_t * endpt,
                              int64_t * naffected,
                              int64_t * affected, int64_t * affected_map,
                              int64_t * ntri, int64_t * delta_global_ntri_in);
extern void sorted_update_tris (int dir, struct stinger *G,
                                int64_t v, int64_t w,
                                int64_t * naffected_in,
                                int64_t * affected,
                                int64_t * ntri,
                                int64_t * delta_global_ntri_in,
                                int64_t * gather_work);
extern void sorted_bulk_update_tris (struct stinger *G, int64_t v,
                                     int64_t nedge, int64_t * endpt,
                                     int64_t * naffected,
                                     int64_t * affected,
                                     int64_t * affected_map, int64_t * ntri,
                                     int64_t * delta_global_ntri_in);
extern void bloom_update_tris (int dir, struct stinger *G, int64_t v,
                               int64_t w, int64_t * naffected_in,
                               int64_t * affected, int64_t * ntri,
                               int64_t * delta_global_ntri_in,
                               uint64_t * filter, int len,
                               int64_t * gather_work);
extern void bloom_bulk_update_tris (struct stinger *G, int64_t v,
                                    int64_t nedge, int64_t * endpt,
                                    int64_t * naffected, int64_t * affected,
                                    int64_t * affected_map, int64_t * ntri,
                                    int64_t * delta_global_ntri_in,
                                    uint64_t * filter, int len);
extern double simple_update_local_cc (struct stinger *G, int64_t naffected,
                                      const int64_t * affected,
                                      const int64_t * ntri, double *local_cc);
extern double simple_update_global_cc (struct stinger *G, int64_t i,
                                       int64_t j, int64_t delta_global_ntri,
                                       double *global_cc,
                                       int64_t * global_ntri,
                                       int64_t * global_degsum);

static int64_t maxdeg, maxdeg2;

static int64_t *restrict affected;

static int64_t global_ntri, global_degsum;
static int64_t global_ntri_init;
static int64_t *restrict ntri;
static int64_t *restrict ntri_init;
static int64_t *restrict ntri_end;

static double global_cc;
static double global_cc_init;
static double *restrict local_cc;
static double *restrict local_cc_init;

static double ntri_time, global_cc_time, local_cc_time, stinger_time,
  update_time;

static int
i64_pair_cmp (const void *pa, const void *pb)
{
  const int64_t ai = ((const int64_t *) pa)[0];
  const int64_t bi = ((const int64_t *) pb)[0];
  if (ai < bi)
    return -1;
  if (ai > bi)
    return 1;
  const int64_t aj = ((const int64_t *) pa)[1];
  const int64_t bj = ((const int64_t *) pb)[1];
  if (aj < bj)
    return -1;
  if (aj > bj)
    return 1;
  return 0;
}


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

  ntri = xcalloc (3 * nv, sizeof (*ntri));
  ntri_end = &ntri[nv];
  ntri_init = &ntri_end[nv];
  local_cc = xcalloc (2 * nv, sizeof (*local_cc));
  local_cc_init = &local_cc[nv];
  affected = xmalloc (nv * sizeof (*affected));

  /* activate threads */
  OMP ("omp parallel for")
    for (int64_t k = 0; k < nv; ++k) {
      affected[k] = 0;
    }

  /* Clustering coefficients */
  OMP ("omp parallel") {
    OMP ("omp master") tic ();
    count_all_triangles (nv, off, ind, ntri);
    OMP ("omp master") ntri_time = toc ();
  }
  PRINT_STAT_DOUBLE ("time_triangles", ntri_time);

  tic ();
  global_ntri = 0;
  OMP ("omp parallel for reduction(+:global_degsum,global_ntri)")
    MTA ("mta assert nodep") MTA ("mta use 100 streams")
    for (int64_t i = 0; i < nv; ++i) {
      const int64_t deg = off[i + 1] - off[i];
      const int64_t d = deg * (deg - 1);
      global_degsum += d;
      global_ntri += ntri[i];
    }
  global_cc = (global_degsum ? global_ntri / (double) global_degsum : 0.0);
  global_cc_time = toc ();
  PRINT_STAT_DOUBLE ("time_global_cc", global_cc_time);

  tic ();
  OMP ("omp parallel for")
    MTA ("mta use 100 streams")
    for (int64_t i = 0; i < nv; ++i) {
      const int64_t deg = off[i + 1] - off[i];
      const int64_t d = deg * (deg - 1);
      local_cc[i] = (d ? ntri[i] / (double) d : 0.0);
    }
  local_cc_time = toc ();
  PRINT_STAT_DOUBLE ("time_local_cc", local_cc_time);

  memcpy (ntri_init, ntri, nv * sizeof (*ntri));
  global_ntri_init = global_ntri;
  global_cc_init = global_cc;
  memcpy (local_cc_init, local_cc, nv * sizeof (*local_cc));

  maxdeg = 0;
  maxdeg2 = 0;
#if !defined(__MTA__)
  OMP ("omp parallel") {
    int64_t lmaxdeg = 0, lmaxdeg2 = 0;
    OMP ("omp for")
      for (int64_t i = 0; i < nv; ++i) {
        const int64_t deg = off[i + 1] - off[i];
        if (deg >= lmaxdeg) {
          lmaxdeg2 = lmaxdeg;
          lmaxdeg = deg;
        } else if (deg > lmaxdeg2) {
          lmaxdeg2 = deg;
        }
      }

    OMP ("omp critical") {
      if (lmaxdeg > maxdeg)
        maxdeg = lmaxdeg;
      if (lmaxdeg2 > maxdeg2)
        maxdeg2 = lmaxdeg2;
    }
  }
#else /* __MTA__ */
  /* The MTA compiler can't handle it, so just use maxdeg twice. */
  for (int64_t i = 0; i < nv; ++i) {
    const int64_t deg = off[i + 1] - off[i];
    if (deg >= maxdeg)
      maxdeg = deg;
  }
  maxdeg2 = maxdeg;
#endif /* __MTA__ */
  /* At most maxdeg + nactions *after* all the actions... */
  maxdeg += nbatch * batch_size;
  if (maxdeg > nv)
    maxdeg = nv;
  maxdeg2 += nbatch * batch_size;
  if (maxdeg2 > nv)
    maxdeg2 = nv;
  /* Round both up to a factor of 16 for cache line happiness. */
  maxdeg = (maxdeg + 15) & ~(int64_t) 15;
  maxdeg2 = (maxdeg2 + 15) & ~(int64_t) 15;


  /* Run the updates */
#if !defined(__MTA__)
  memcpy (ntri, ntri_init, nv * sizeof (*ntri));
  memcpy (local_cc, local_cc_init, nv * sizeof (*local_cc));
#else
  MTA ("mta use 100 streams")
    for (size_t k = 0; k < nv; ++k) {
      ntri[k] = ntri_init[k];
    }

  MTA ("mta use 100 streams")
    for (size_t k = 0; k < nv; ++k) {
      local_cc[k] = local_cc_init[k];
    }
#endif
  global_ntri = global_ntri_init;
  global_cc = global_cc_init;

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  stinger_time = toc ();
  PRINT_STAT_DOUBLE ("time_stinger", stinger_time);

  int ndis = 0;
  int64_t *off_end = NULL, *ind_end = NULL, *ewgt_end = NULL;

  stinger_to_sorted_csr (S, nv, &off_end, &ind_end, &ewgt_end, NULL, NULL, NULL);
  count_all_triangles (nv, off_end, ind_end, &ntri_end[0]);

  OMP ("omp parallel for reduction(+:ndis)")
    MTA ("mta assert parallel") MTA ("mta use 100 streams")
    for (int64_t i = 0; i < nv; ++i) {
      if (ntri[i] != ntri_end[i] && ndis < MAX_DISAGREE_PRINTED) {
        ++ndis;
        fprintf (stderr, "ERROR: vertex %ld disagrees, %ld != %ld\n",
                 (long int) i, (long int) ntri[i], (long int) ntri_end[i]);
      }
    }
  free (ewgt_end);
  free (ind_end);
  free (off_end);
  free (graphmem);

  tic ();

  int64_t ntrace = 0;
  if (1 == batch_size) {
    int64_t gather_work[maxdeg + maxdeg2];

    for (int actk = 0; actk < nbatch; ++actk, ++ntrace) {
      double t;
      double global_delta, mld = 0.0, gd = 0.0;
      int64_t delta_global_ntri = 0, naffected = 0;
      int changed;

      int64_t i = ACTI (actk);
      int64_t j = ACTJ (actk);

      if (i < 0) {
        /* Delete. */
        i = -i - 1;
        j = -j - 1;

        assert (i < nv);
        assert (j < nv);
        assert (i >= 0);
        assert (j >= 0);

        changed = stinger_remove_edge_pair (S, 0, i, j);

        if (changed && i != j) {
          sorted_update_tris (-1, S, i, j, &naffected, &affected[0],
                              &ntri[0], &delta_global_ntri, gather_work);
          mld = simple_update_local_cc (S, naffected, affected, ntri,
                                        &local_cc[0]);
          gd = simple_update_global_cc (S, i, j, delta_global_ntri,
                                        &global_cc,
                                        &global_ntri, &global_degsum);
        }
      } else {                    /* Add an edge. */
        assert (i < nv);
        assert (j < nv);
        assert (i >= 0);
        assert (j >= 0);

        changed = stinger_insert_edge_pair (S, 0, i, j, 1, ntrace + 1);

        if (changed && i != j) {
          sorted_update_tris (1, S, i, j, &naffected, &affected[0],
                              &ntri[0], &delta_global_ntri, gather_work);
          mld = simple_update_local_cc (S, naffected, affected, ntri,
                                        &local_cc[0]);
          gd = simple_update_global_cc (S, i, j, delta_global_ntri,
                                        &global_cc,
                                        &global_ntri, &global_degsum);
        }
      }
    }
  } else {
    double t;
    int64_t naffected = 0;
    int64_t * affected_map = xmalloc (nv * sizeof(*affected_map));

    for (int64_t k = 0; k < nv; ++k) {
      affected_map[k] = -1;
    }

    int64_t * act = xmalloc (2 * 2 * batch_size * sizeof(*act));
    int * insoff = xmalloc ( (2 * 2 * batch_size + 1) * sizeof(*insoff));
    int * deloff = xmalloc ( (2 * 2 * batch_size + 1) * sizeof(*deloff));
    int * actk = xmalloc ( (2 * 2 * batch_size + 1) * sizeof(*actk));

    for (int actno = 0; actno < nbatch * batch_size; actno += batch_size) {
      const int endact = (actno + batch_size > naction ?
                          naction : actno + batch_size);
      int nact = 2 * (endact - actno);
      double mld = 0.0, gd = 0.0;
      int64_t delta_global_ntri = 0;

      /*
        If a deletion comes before an insertion, and the
        edge is *not* in the graph, this will be wrong.
        Our test cases always choose deletions from
        generated edges, so we can ignore the order for
        now.

        Sort the actions by I and encoded J.
        Associate each I with an index.
        Compute add, del offsets.
        Convert del indices.
      */
      int n = 1;
      int k2 = 0;

      OMP ("omp parallel") {
        /* Copy & make i positive if necessary.  Negative j
           still implies deletion. */
        OMP ("omp for")
          MTA ("mta assert nodep")
          MTA ("mta block schedule")
          for (int k = actno; k < endact; ++k) {
            const int64_t i = ACTI (k);
            const int64_t j = ACTJ (k);
            if (i != j) {         /* Skip self edges */
              int where = stinger_int_fetch_add (&k2, 2);
              assert (where < 2 * batch_size);
              ACTI2 (where) = (i < 0 ? -i - 1 : i);
              ACTJ2 (where) = j;
              ACTI2 (where + 1) = (j < 0 ? -j - 1 : j);
              ACTJ2 (where + 1) = i;
            }
          }

        OMP ("omp single") {
          nact = k2;
          qsort (&ACTI2 (0), nact, 2 * sizeof (act[0]), i64_pair_cmp);
          actk[0] = 0;
        }

        /* Find local indices... */
        OMP ("omp for") MTA ("mta assert nodep")
          for (int k = 1; k < nact; ++k) {
            if (ACTI2 (k) != ACTI2 (k - 1)) {
              stinger_int_fetch_add (&n, 1);
              actk[k] = 1;
            } else {
              actk[k] = 0;
            }
          }

        OMP ("omp single")      // yes, I need to implement this...
          for (int k = 1; k < nact; ++k) {
            actk[k] += actk[k - 1];
          }
        assert (actk[nact - 1] == n - 1);

        /* Offsets. */
        OMP ("omp for") MTA ("mta block schedule")
          for (int k = 0; k <= n; ++k) {
            deloff[k] = 0;
          }

        OMP ("omp for") MTA ("mta block schedule")
          for (int k = 0; k < nact; ++k) {
            stinger_int_fetch_add (&deloff[actk[k] + 1], 1);
          }

        OMP ("omp single")
          for (int k = 1; k <= n; ++k) {
            deloff[k] += deloff[k - 1];
          }
        assert (deloff[n] == nact);

        OMP ("omp for") MTA ("mta assert nodep")
          MTA ("mta loop norestructure")
          for (int k = 0; k < n; ++k) {
            int off;
            const int endoff = deloff[k + 1];

            for (off = deloff[k]; off < endoff; ++off) {
              if (ACTJ2 (off) >= 0)
                break;
            }
            insoff[k] = off;
            assert (insoff[k] == deloff[k + 1]
                    || ACTI2 (insoff[k]) == ACTI2 (deloff[k]));

            for (off = deloff[k]; off < insoff[k]; ++off) {
              ACTJ2 (off) = -ACTJ2 (off) - 1;
            }
          }

        OMP ("omp for")
	  for (int k = 0; k < n; ++k) {
            const int64_t i = ACTI2 (deloff[k]);
            const int64_t deg = stinger_outdegree (S, i);
            stinger_int64_fetch_add (&global_degsum, -deg * (deg - 1));
          }

        OMP ("omp for") MTA ("mta assert nodep") MTASTREAMS ()
          for (int k = 0; k < n; ++k) {
            const int64_t i = ACTI2 (deloff[k]);
            const int64_t nrem = insoff[k] - deloff[k];
            const int64_t nins = deloff[k + 1] - insoff[k];

            int64_t *restrict torem, *restrict toins;
            torem =
              (nrem + nins ? malloc ((nrem + nins) * sizeof (*torem)) : NULL);
            toins = (torem ? &torem[nrem] : NULL);

            MTA ("mta assert nodep") MTA ("mta block schedule")
              for (int k2 = 0; k2 < nrem; ++k2) {
                torem[k2] = ACTJ2 (deloff[k] + k2);
                assert (ACTI2 (deloff[k] + k2) == i);
              }

            MTA ("mta assert nodep") MTA ("mta block schedule")
              for (int k2 = 0; k2 < nins; ++k2) {
                toins[k2] = ACTJ2 (insoff[k] + k2);
                assert (ACTI2 (deloff[k] + k2) == i);
              }

            stinger_remove_and_insert_edges (S, 0, i,
                                             nrem, torem, nins, toins,
                                             NULL, actno);

            free (torem);
          }

        OMP ("omp for")
	  for (int k = 0; k < n; ++k) {
            const int64_t i = ACTI2 (deloff[k]);
            const int64_t deg = stinger_outdegree (S, i);
            stinger_int64_fetch_add (&global_degsum, +deg * (deg - 1));
          }

        /* Assume each changed. */
        /* Scatter/gather affected lists */
        OMP ("omp for")
          MTA ("mta assert nodep") MTA ("mta block dynamic schedule")
          MTASTREAMS ()
          for (int k = 0; k < n; ++k) {
            const int nvtx = deloff[k + 1] - deloff[k];
            int64_t vtx[nvtx];
            const int64_t i = ACTI2 (deloff[k]);
            int vk = 0, k2;

            for (k2 = deloff[k]; k2 < insoff[k]; ++k2, ++vk) {
              vtx[vk] = ~ACTJ2 (k2);
            }

            for (; k2 < deloff[k + 1]; ++k2, ++vk) {
              vtx[vk] = ACTJ2 (k2);
            }

            sorted_bulk_update_tris (S, i, nvtx, vtx,
                                     &naffected, affected, affected_map,
                                     &ntri[0], &delta_global_ntri);
          }

        mld = simple_update_local_cc (S, naffected, affected, ntri,
                                      &local_cc[0]);

        OMP ("omp single") {
          global_ntri += delta_global_ntri;
          double new_global_cc =
            (global_degsum ? global_ntri / (double) global_degsum : 0.0);
          gd = global_cc - new_global_cc;
          global_cc = new_global_cc;
        }

        OMP ("omp master") {
          ++ntrace;
        }

        OMP ("omp for") MTA ("mta assert nodep")
          for (int64_t k = 0; k < naffected; ++k) {
            affected_map[affected[k]] = -1;
          }
        naffected = 0;
        OMP ("omp barrier")}
    }

    free (actk);
    free (insoff);
    free (deloff);
    free (act);
    free (affected_map);
  }

  update_time = toc ();
  PRINT_STAT_DOUBLE ("time_updates", update_time);

  free (affected);
  free (ntri);
  free (local_cc);
  stinger_free_all (S);
  free (actionmem);

  printf ("\n}\n");
}

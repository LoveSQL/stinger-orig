#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include <assert.h>

#include "stinger.h"
#include "stinger-atomics.h"

#if defined(__MTA__)
#define MTA(x) _Pragma(x)
#if defined(MTA_STREAMS)
#define MTASTREAMS() MTA(MTA_STREAMS)
#else
#define MTASTREAMS() MTA("mta use 100 streams")
#endif
#else
#define MTA(x)
#define MTASTREAMS()
#endif

#if defined(_OPENMP)
#define OMP(x) _Pragma(x)
#else
#define OMP(x)
#endif

MTA ("mta inline")
static inline void
add_to_affected (const int64_t v,
                 int64_t * restrict naffected,
                 int64_t * restrict affected,
                 int64_t * restrict affected_map)
{
  if (affected_map[v] >= 0)
    return;

#if !defined(__MTA__)
  if (-1 != stinger_int64_cas (&affected_map[v], -1, INT64_MAX))
    /* Already claimed by someone else. */
    return;

  int64_t which = stinger_int64_fetch_add (naffected, 1);
  affected[which] = v;
  affected_map[v] = which;
#else
  int64_t where = readfe (&affected_map[v]);
  if (where < 0) {
    where = int_fetch_add (naffected, 1);
    affected[where] = v;
  }
  writeef (&affected_map[v], where);
#endif
}

void
simple_update_tris (int dir, struct stinger *G, int64_t v, int64_t w,
                    size_t * restrict naffected_in,
                    int64_t * restrict affected, int64_t * restrict ntri,
                    int64_t * restrict delta_global_ntri_in,
                    int64_t * restrict gather_work)
{
  size_t vdeg, vdegmax, wdeg, wdegmax;
  size_t *vdegptr = &vdeg;
  size_t *wdegptr = &wdeg;
  int64_t *restrict vn = NULL;
  int64_t *restrict wn = NULL;
  size_t naffected = 2;
  int64_t delta_global_ntri = 0;
  int64_t incr_v = 0, incr_w = 0;
  int64_t incr = 0;

  if (v == w) {
    if (affected)
      *naffected_in = 0;
    return;
  }

  if (dir < 0)
    dir = -1;
  else
    dir = 1;

  /* By this point, they at least have their degrees affected. */
  if (affected)
    affected[0] = v;
  if (affected)
    affected[1] = w;
  /* Gather neighbor lists */
  vdegmax = stinger_outdegree (G, v);
  wdegmax = stinger_outdegree (G, w);
  if (!wdegmax || !vdegmax) {
    *delta_global_ntri_in = 0;
    if (affected)
      *naffected_in = 0;
    return;
  }
  vn = gather_work;
  vdegmax = (vdegmax + 15) & ~(size_t) 15;      /* Rounded for cache-line foo. */
  wn = &vn[vdegmax];
  OMP ("omp parallel") {
#if defined(__MTA__)
    {
      volatile int tmp;
      future int v$, w$;
      future v$ (G, v, vdegptr, vn, vdegmax) {
        stinger_gather_typed_successors (G, 0, v, vdegptr, vn, vdegmax);
        return 1;
      }
      future w$ (G, w, wdegptr, wn, wdegmax) {
        stinger_gather_typed_successors (G, 0, w, wdegptr, wn, wdegmax);
        return 1;
      }
      tmp = v$ + w$;
    }
#else
    OMP ("omp sections") {
      OMP ("omp section") stinger_gather_typed_successors (G, 0, v, &vdeg, vn,
                                                           vdegmax);
      OMP ("omp section") stinger_gather_typed_successors (G, 0, w, &wdeg, wn,
                                                           wdegmax);
    }
#endif

    assert (vdeg == stinger_outdegree (G, v));
    assert (wdeg == stinger_outdegree (G, w));

    /* Karl's idea: scan larger for smaller */
    if (wdeg && vdeg >= wdeg) {
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        MTA ("mta assert nodep")
        for (size_t kv = 0; kv < vdeg; ++kv) {
          const int64_t neigh = vn[kv];
          if (neigh != v)
            MTA ("mta assert nodep")
              for (size_t kw = 0; kw < wdeg; ++kw) {
                if (w != wn[kw] && neigh == wn[kw]) {
                  incr += dir;
                  stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                  if (affected) {
                    const size_t kaffected =
                      stinger_size_fetch_add (&naffected, 1);
                    affected[kaffected] = neigh;
                  }
#if !defined(__MTA__)
                  break;
#endif
                }
              }
        }
    } else if (vdeg) {
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        for (size_t kw = 0; kw < wdeg; ++kw) {
          const int64_t neigh = wn[kw];
          if (neigh != w)
            MTA ("mta assert nodep")
              for (size_t kv = 0; kv < vdeg; ++kv) {
                if (v != vn[kv] && neigh == vn[kv]) {
                  incr += dir;
                  stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                  if (affected) {
                    const size_t kaffected =
                      stinger_size_fetch_add (&naffected, 1);
                    affected[kaffected] = neigh;
                  }
#if !defined(__MTA__)
                  break;
#endif
                }
              }
        }
    }
  }

  ntri[v] += 2 * incr;
  ntri[w] += 2 * incr;
  *delta_global_ntri_in = 6 * incr;
  if (affected)
    *naffected_in = naffected;
}

void
bulk_update_tris (struct stinger *G, int64_t v,
                  int64_t nedge, int64_t * restrict endpt,
                  int64_t * restrict naffected,
                  int64_t * restrict affected,
                  int64_t * restrict affected_map, int64_t * restrict ntri,
                  int64_t * restrict delta_global_ntri_in)
{
  size_t vdeg, vdegmax;
  int64_t delta_global_ntri = 0;
  int64_t incr_v = 0, incr_w = 0;
  int64_t incr = 0;

  if (affected)
    add_to_affected (v, naffected, affected, affected_map);
  vdegmax = stinger_outdegree (G, v);
  if (!vdegmax) {
    *delta_global_ntri_in = 0;
    if (affected)
      *naffected = 0;
    return;
  }
  int64_t vn[vdegmax];
  stinger_gather_typed_successors (G, 0, v, &vdeg, vn, vdegmax);
#if 1
  if (vdeg != stinger_outdegree (G, v))
    fprintf (stderr, "%d: vdeg %d  outdeg %d\n", (int) v, (int) vdeg,
             (int) stinger_outdegree (G, v));
#endif
  assert (vdeg == stinger_outdegree (G, v));

  MTA ("mta assert parallel") MTASTREAMS ()
    for (int64_t k = 0; k < nedge; ++k) {
      size_t wdeg;
      const int dir = (endpt[k] < 0 ? -1 : 1);
      const int64_t w = (endpt[k] < 0 ? ~endpt[k] : endpt[k]);
      const int64_t wdegmax = stinger_outdegree (G, w);
      //int64_t wn[wdegmax];
      int64_t *restrict wn;
      if (!wdegmax)
        continue;
      wn = malloc (wdegmax * sizeof (*wn));
      stinger_gather_typed_successors (G, 0, w, &wdeg, wn, wdegmax);
      if (affected)
        add_to_affected (w, naffected, affected, affected_map);

      /* Karl's idea: scan larger for smaller */
      if (wdeg && vdeg >= wdeg) {
        //OMP("omp for reduction(+:incr)")
        MTA ("mta loop restructure")
          MTA ("mta assert nodep")
          for (size_t kv = 0; kv < vdeg; ++kv) {
            const int64_t neigh = vn[kv];
            if (neigh != v)
              MTA ("mta assert nodep")
                for (size_t kw = 0; kw < wdeg; ++kw) {
                  if (w != wn[kw] && neigh == wn[kw]) {
                    incr += dir;
                    stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                    if (affected)
                      add_to_affected (neigh, naffected, affected, affected_map);
#if !defined(__MTA__)
                    break;
#endif
                  }
                }
          }
      } else if (vdeg) {
        //OMP("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
          for (size_t kw = 0; kw < wdeg; ++kw) {
            const int64_t neigh = wn[kw];
            if (neigh != w)
              MTA ("mta assert nodep")
                for (size_t kv = 0; kv < vdeg; ++kv) {
                  if (v != vn[kv] && neigh == vn[kv]) {
                    incr += dir;
                    stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                    if (affected)
                      add_to_affected (neigh, naffected, affected, affected_map);
#if !defined(__MTA__)
                    break;
#endif
                  }
                }
          }
      }
      free (wn);
      stinger_int64_fetch_add (&ntri[v], 2 * incr);
      stinger_int64_fetch_add (&ntri[w], 2 * incr);
      stinger_int64_fetch_add (&delta_global_ntri, 6 * incr);
    }

  stinger_int64_fetch_add (delta_global_ntri_in, delta_global_ntri);
}

double
simple_update_local_cc (struct stinger *G,
                        size_t naffected, const int64_t * restrict affected,
                        const int64_t * restrict ntri,
                        double *restrict local_cc)
{
  double maxmagdiff = 0.0, outdiff = 0.0;
  OMP ("omp parallel") {
    double lmax = 0.0;
    MTA ("mta assert nodep *local_cc")
      OMP ("omp for")
      for (size_t k = 0; k < naffected; ++k) {
        const int64_t v = affected[k];
        const double d = stinger_outdegree (G, v);
        const double new_local_cc = (d > 1 ? ntri[v] / (d * (d - 1)) : 0.0);
        const double diff = new_local_cc - local_cc[v];
        const double diffmag = fabs (diff);
        local_cc[v] = new_local_cc;
        if (diffmag > lmax) {
          lmax = diffmag;
        }
        //if (diffmag > maxmagdiff) { outdiff = diff; }
      }
    if (lmax > maxmagdiff) {
      OMP ("omp critical") {
        if (lmax > maxmagdiff)
          maxmagdiff = lmax;
      }
    }
  }
  //return outdiff;
  return maxmagdiff;
}

double
simple_update_global_cc (struct stinger *G,
                         int64_t i, int64_t j,
                         int64_t delta_global_ntri,
                         double *restrict global_cc,
                         int64_t * restrict global_ntri,
                         int64_t * restrict global_degsum)
{
  /* These are the degrees *after* adding/removing the edge. */
  const int64_t ideg = stinger_outdegree (G, i);
  const int64_t jdeg = stinger_outdegree (G, j);
  int64_t gntri, gdegsum, incr;
  double new_global_cc, out;

  gntri = stinger_int64_fetch_add (global_ntri, delta_global_ntri)
    + delta_global_ntri;
  incr = 2 * (ideg + jdeg - 2);
  if (delta_global_ntri < 0)    /* removed an edge */
    incr = -incr;
  gdegsum = stinger_int64_fetch_add (global_degsum, incr) + incr;
  new_global_cc = (gdegsum ? gntri / (double) gdegsum : 0.0);
  out = new_global_cc - *global_cc;
  *global_cc = new_global_cc;
  return out;
}

static int
i64cmp (const void *ap, const void *bp)
{
  return (int) ((*(int64_t *) ap) - (*(int64_t *) bp));
}

extern size_t count_intersections (const int64_t ai, const size_t alen,
                                   const int64_t * restrict a,
                                   const int64_t bi, const size_t blen,
                                   const int64_t * restrict b);

#define SEQLEN 8

MTA ("mta inline") MTA ("mta expect parallel context")
static
int64_t
find_in_sorted (const int64_t tofind,
                const int64_t N, const int64_t * restrict ary)
{
  int64_t out = -1;
  if (N <= 0)
    return -1;

  int64_t top = N - 1, bot = 0;

  if (tofind == ary[bot])
    return bot;
  if (tofind == ary[top])
    return top;
  while (top - bot + 1 > SEQLEN) {
    int64_t mid = bot + ((top - bot) / 2);
    if (tofind == ary[mid])
      return mid;
    if (tofind < ary[mid])
      bot = mid;
    else
      top = mid;
  }
  for (; bot <= top; ++bot)
    if (tofind == ary[bot])
      return bot;
  return -1;
}

void
sorted_update_tris (int dir, struct stinger *G, int64_t v, int64_t w,
                    size_t * restrict naffected_in,
                    int64_t * restrict affected, int64_t * restrict ntri,
                    int64_t * restrict delta_global_ntri_in,
                    int64_t * restrict gather_work)
{
  size_t vdeg, vdegmax, wdeg, wdegmax;
  size_t *vdegptr = &vdeg;
  size_t *wdegptr = &wdeg;
  int64_t *restrict vn = NULL;
  int64_t *restrict wn = NULL;
  size_t naffected = 2;
  int64_t delta_global_ntri = 0;
  int64_t incr_v = 0, incr_w = 0;
  int64_t incr = 0;

  if (v == w) {
    if (affected)
      *naffected_in = 0;
    return;
  }

  if (dir < 0)
    dir = -1;
  else
    dir = 1;

  /* By this point, they at least have their degrees affected. */
  if (affected)
    affected[0] = v;
  if (affected)
    affected[1] = w;
  /* Gather neighbor lists */
  vdegmax = stinger_outdegree (G, v);
  wdegmax = stinger_outdegree (G, w);
  if (!wdegmax || !vdegmax) {
    *delta_global_ntri_in = 0;
    if (affected)
      *naffected_in = 0;
    return;
  }
  vn = gather_work;
  wn = &vn[(vdegmax + 15) & ~(size_t) 15];      /* Rounded for cache-line foo. */
  OMP ("omp parallel") {
#if defined(__MTA__)
    {
      volatile int tmp;
      future int v$, w$;
      future v$ (G, v, vdegptr, vn, vdegmax) {
        stinger_gather_typed_successors (G, 0, v, vdegptr, vn, vdegmax);
        return 1;
      }
      future w$ (G, w, wdegptr, wn, wdegmax) {
        stinger_gather_typed_successors (G, 0, w, wdegptr, wn, wdegmax);
        return 1;
      }
      tmp = v$ + w$;
    }
#else
    OMP ("omp sections") {
      OMP ("omp section") stinger_gather_typed_successors (G, 0, v, &vdeg, vn,
                                                           vdegmax);
      OMP ("omp section") stinger_gather_typed_successors (G, 0, w, &wdeg, wn,
                                                           wdegmax);
    }
#endif

    assert (vdeg == stinger_outdegree (G, v));
    assert (wdeg == stinger_outdegree (G, w));

    if (vdeg && vdeg <= wdeg) {
      OMP ("omp single") qsort (&wn[0], wdeg, sizeof (wn[0]), i64cmp);
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        for (size_t kv = 0; kv < vdeg; ++kv) {
          const int64_t neigh = vn[kv];
          if (neigh != v && neigh != w) {
            /* Binary search while large. */
            size_t high = wdeg - 1, low = 0;
            int found = 0;

            while (high > low + 8) {
              const size_t pt = (high + low) >> 1;
              if (neigh == wn[pt]) {
                incr += dir;
                stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                if (affected) {
                  const size_t kaffected =
                    stinger_size_fetch_add (&naffected, 1);
                  affected[kaffected] = neigh;
                }
                found = 1;
                break;
              } else if (neigh > wn[pt])
                low = pt;
              else
                high = pt;
            }
            if (!found)
              MTA ("mta assert nodep")
                MTA ("mta interleave schedule")
                for (size_t k = low; k <= high; ++k)
                  if (neigh == wn[k]) {
                    incr += dir;
                    stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                    if (affected) {
                      const size_t kaffected =
                        stinger_size_fetch_add (&naffected, 1);
                      affected[kaffected] = neigh;
                    }
#if !defined(__MTA__)
                    break;
#endif
                  }
          }
        }
    } else if (wdeg) {
      OMP ("omp single") qsort (&vn[0], vdeg, sizeof (vn[0]), i64cmp);
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        for (size_t kw = 0; kw < wdeg; ++kw) {
          const int64_t neigh = wn[kw];
          if (neigh != w && neigh != v) {
            /* Binary search. */
            size_t high = vdeg - 1, low = 0;
            int found = 0;

            while (high > low + 8) {
              const size_t pt = (high + low) >> 1;
              if (neigh == vn[pt]) {
                incr += dir;
                stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                if (affected) {
                  const size_t kaffected =
                    stinger_size_fetch_add (&naffected, 1);
                  affected[kaffected] = neigh;
                }
                found = 1;
                break;
              } else if (neigh > vn[pt])
                low = pt;
              else
                high = pt;
            }
            if (!found)
              MTA ("mta assert nodep")
                MTA ("mta interleave schedule")
                for (size_t k = low; k <= high; ++k)
                  if (neigh == vn[k]) {
                    incr += dir;
                    stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
                    if (affected) {
                      const size_t kaffected =
                        stinger_size_fetch_add (&naffected, 1);
                      affected[kaffected] = neigh;
                    }
#if !defined(__MTA__)
                    break;
#endif
                  }
          }
        }
    }
  }

  ntri[v] += 2 * incr;
  ntri[w] += 2 * incr;
  *delta_global_ntri_in = 6 * incr;
  if (affected)
    *naffected_in = naffected;
}

void
sorted_bulk_update_tris (struct stinger *G, int64_t v,
                         int64_t nedge, int64_t * restrict endpt,
                         int64_t * restrict naffected,
                         int64_t * restrict affected,
                         int64_t * restrict affected_map,
                         int64_t * restrict ntri,
                         int64_t * restrict delta_global_ntri_in)
{
  size_t vdeg, vdegmax;
  int64_t delta_global_ntri = 0;

  if (affected)
    add_to_affected (v, naffected, affected, affected_map);
  vdegmax = stinger_outdegree (G, v);
  if (!vdegmax) {
    *delta_global_ntri_in = 0;
    if (affected)
      *naffected = 0;
    return;
  }
  int64_t vn[vdegmax];
  stinger_gather_typed_successors (G, 0, v, &vdeg, vn, vdegmax);
  assert (vdeg == stinger_outdegree (G, v));
  qsort (&vn[0], vdeg, sizeof (vn[0]), i64cmp);

  MTA ("mta assert parallel") MTASTREAMS ()
    for (int64_t k = 0; k < nedge; ++k) {
      int64_t incr = 0;
      const int dir = (endpt[k] < 0 ? -1 : 1);
      const int64_t w = (endpt[k] < 0 ? ~endpt[k] : endpt[k]);

      if (k > 0 && endpt[k] == endpt[k-1]) continue;

      const int64_t wdegmax = stinger_outdegree (G, w);
      size_t wdeg;
      //int64_t wn[wdegmax];
      int64_t *restrict wn;
      if (!wdegmax)
        continue;
      wn = malloc (wdegmax * sizeof (*wn));
      stinger_gather_typed_successors (G, 0, w, &wdeg, wn, wdegmax);
      assert (wdeg == stinger_outdegree (G, w));
      if (affected)
        add_to_affected (w, naffected, affected, affected_map);

      /* Search. */
      MTA ("mta assert nodep")
        for (int64_t k2 = 0; k2 < wdeg; ++k2) {
          const int64_t where = find_in_sorted (wn[k2], vdeg, vn);
          if (where >= 0) {
            incr += dir;
            stinger_int64_fetch_add (&ntri[wn[k2]], 2 * dir);
            if (affected)
              add_to_affected (wn[k2], naffected, affected, affected_map);
          }
        }

      stinger_int64_fetch_add (&ntri[v], 2 * incr);
      stinger_int64_fetch_add (&ntri[w], 2 * incr);
      stinger_int64_fetch_add (&delta_global_ntri, 6 * incr);
      free (wn);
    }

  stinger_int64_fetch_add (delta_global_ntri_in, delta_global_ntri);
}

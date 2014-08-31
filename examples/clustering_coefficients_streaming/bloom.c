
#include <limits.h>
#include "stinger.h"
#include "stinger-atomics.h"


MTA ("mta inline") MTA ("mta expect parallel context")
static
uint64_t
hash64shift (uint64_t key)
{
  key = (~key) + (key << 21);   // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8);        // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4);        // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

MTA ("mta inline") MTA ("mta expect parallel context")
static int
bloom_insert (uint64_t * filter, int len, int key)
{
  const unsigned hash = hash64shift (key) % (len << 6);
  const unsigned word_offset = hash >> 6;
  const unsigned bit_offset = hash & 0x3Fu;
  const uint64_t bit = ((uint64_t) 1) << bit_offset;

  //assert (word_offset < len);

  OMP ("omp atomic") MTA ("mta update") filter[word_offset] |= bit;
  //printf("key: %d, hash: %d, offset: %d, bit: %d\n", key, hash, word_offset, bit_offset);
  return 0;
}

MTA ("mta inline") MTA ("mta expect parallel context")
static
_Bool
bloom_query (uint64_t * filter, int len, int key)
{
  const unsigned hash = hash64shift (key) % (len << 6);
  const unsigned word_offset = hash >> 6;
  const unsigned bit_offset = hash & 0x3Fu;
  const uint64_t bit = ((uint64_t) 1) << bit_offset;

  //assert (word_offset < len);

  return bit & filter[word_offset];

#if 0
  uint64_t test;
  if (bit_offset == 63)
    test = 0x7FFFFFFFFFFFFFFFull;
  else {
    test = 0xBFFFFFFFFFFFFFFFull;
    test = test >> (62 - bit_offset);
  }

  //    printf("\nkey: %d, hash: %d, offset: %d, bit: %d\n", key, hash, word_offset, bit_offset);
  //    printf("test: %lld\n", test);
  if (~(test | filter[word_offset]) == 0)
    return 1;
  else
    return 0;
#endif
}

void
bloom_update_tris (int dir, struct stinger *G, int64_t v, int64_t w,
                   int64_t * restrict naffected_in,
                   int64_t * restrict affected, int64_t * restrict ntri,
                   int64_t * restrict delta_global_ntri_in,
                   uint64_t * restrict filter, int len,
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

  assert (sizeof (uint64_t) * CHAR_BIT == 64);

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
      OMP ("omp section")
        stinger_gather_typed_successors (G, 0, v, &vdeg, vn, vdegmax);
      OMP ("omp section")
        stinger_gather_typed_successors (G, 0, w, &wdeg, wn, wdegmax);
    }
#endif

    assert (vdeg == stinger_outdegree (G, v));
    assert (wdeg == stinger_outdegree (G, w));

    OMP ("omp for")
      for (int k = 0; k < len; ++k)
        filter[k] = 0;

    /* Karl's idea: scan larger for smaller */
    if (vdeg <= wdeg) {
      OMP ("omp for")
        MTA ("mta assert nodep")
        for (int i = 0; i < vdeg; i++)
          if (vn[i] != v)
            bloom_insert (filter, len, vn[i]);
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        for (size_t kw = 0; kw < wdeg; ++kw) {
          if (w != wn[kw] && bloom_query (filter, len, wn[kw])) {
            const int64_t neigh = wn[kw];
            incr += dir;
            stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
            if (affected) {
              const size_t kaffected = stinger_size_fetch_add (&naffected, 1);
              affected[kaffected] = neigh;
            }
#if !defined(__MTA__)&&!defined(_OPENMP)
            break;
#endif
          }
        }
    } else {
      OMP ("omp for")
        MTA ("mta assert nodep")
        for (int i = 0; i < wdeg; i++)
          if (wn[i] != w)
            bloom_insert (filter, len, wn[i]);
      OMP ("omp for reduction(+:incr)")
        MTA ("mta assert nodep")
        for (size_t kv = 0; kv < vdeg; ++kv) {
          if (v != vn[kv] && bloom_query (filter, len, vn[kv])) {
            const int64_t neigh = vn[kv];
            incr += dir;
            stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
            if (affected) {
              const size_t kaffected = stinger_size_fetch_add (&naffected, 1);
              affected[kaffected] = neigh;
            }
#if !defined(__MTA__)&&!defined(_OPENMP)
            break;
#endif
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
bloom_bulk_update_tris (struct stinger *G, int64_t v,
                        int64_t nedge, int64_t * restrict endpt,
                        int64_t * restrict naffected,
                        int64_t * restrict affected,
                        int64_t * restrict affected_map,
                        int64_t * restrict ntri,
                        int64_t * restrict delta_global_ntri_in,
                        uint64_t * restrict filter, int len)
{
  size_t vdeg, vdegmax;
  int64_t delta_global_ntri = 0;

  /* By this point, they at least have their degrees affected. */
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
  for (int k = 0; k < len; ++k)
    filter[k] = 0;
  for (int i = 0; i < vdeg; i++)
    if (vn[i] != v)
      bloom_insert (filter, len, vn[i]);

  MTA ("mta assert parallel")
    for (int64_t k = 0; k < nedge; ++k) {
      int64_t incr = 0;
      const int dir = (endpt[k] < 0 ? -1 : 1);
      const int64_t w = (endpt[k] < 0 ? ~endpt[k] : endpt[k]);
      const int64_t wdegmax = stinger_outdegree (G, w);
      size_t wdeg;
      //int64_t wn[wdegmax];
      int64_t *restrict wn;
      if (!wdegmax)
        continue;
      wn = malloc (wdegmax * sizeof (*wn));
      stinger_gather_typed_successors (G, 0, w, &wdeg, wn, wdegmax);
      if (affected)
        add_to_affected (w, naffected, affected, affected_map);

      /* Search. */
      MTA ("mta assert nodep")
        for (size_t kw = 0; kw < wdeg; ++kw) {
          const int64_t neigh = wn[kw];
          if (w != neigh && bloom_query (filter, len, neigh)) {
            incr += dir;
            stinger_int64_fetch_add (&ntri[neigh], 2 * dir);
            if (affected)
              add_to_affected (neigh, naffected, affected, affected_map);
          }
        }

      stinger_int64_fetch_add (&ntri[v], 2 * incr);
      stinger_int64_fetch_add (&ntri[w], 2 * incr);
      stinger_int64_fetch_add (&delta_global_ntri, 6 * incr);
      free (wn);
    }

  stinger_int64_fetch_add (delta_global_ntri_in, delta_global_ntri);
}

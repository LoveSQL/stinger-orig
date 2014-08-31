#include <stdlib.h>
#include <inttypes.h>

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

#if defined(__GNUC__)
#define UNLIKELY(x) __builtin_expect ((x), 0)
#else
#define UNLIKELY(x) (x)
#endif

static inline int64_t count_triangles (const size_t nv,
                                       const int64_t * restrict off,
                                       const int64_t * restrict ind,
                                       int64_t i);
static inline size_t count_intersections (const int64_t ai, const size_t alen,
                                          const int64_t * restrict a,
                                          const int64_t bi, const size_t blen,
                                          const int64_t * restrict b);

MTA ("mta expect parallel context") MTA ("mta inline") size_t
count_intersections (const int64_t ai, const size_t alen,
                     const int64_t * restrict a, const int64_t bi,
                     const size_t blen, const int64_t * restrict b)
{
  size_t ka = 0, kb = 0;
  size_t out = 0;

  if (!alen || !blen || a[alen - 1] < b[0] || b[blen - 1] < a[0])
    return 0;

  while (1) {
    if (ka >= alen || kb >= blen)
      break;

    int64_t va = a[ka];
    int64_t vb = b[kb];

    /* Skip self-edges. */
    if (UNLIKELY (va == ai)) {
      ++ka;
      if (ka >= alen)
        break;
      va = a[ka];
    }
    if (UNLIKELY (vb == bi)) {
      ++kb;
      if (kb >= blen)
        break;
      vb = b[kb];
    }

    if (va == vb) {
      ++ka;
      ++kb;
      ++out;
    } else if (va < vb) {
      /* Advance ka */
      ++ka;
      while (ka < alen && a[ka] < vb)
        ++ka;
    } else {
      /* Advance kb */
      ++kb;
      while (kb < blen && va > b[kb])
        ++kb;
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
  //OMP("omp parallel firstprivate(off, ind, ntri)")
  {
    const size_t N = nv;
    OMP ("omp for schedule(dynamic,128)")
      MTA ("mta dynamic schedule") MTASTREAMS ()
      for (size_t i = 0; i < N; ++i)
        ntri[i] = count_triangles (nv, off, ind, i);
  }
}

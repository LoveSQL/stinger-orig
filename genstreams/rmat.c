/* -*- mode: C; mode: folding; fill-column: 70; -*- */

#include "stinger.h"
#include "stinger-atomics.h"
#include "xmalloc.h"
#include "prng.h"


static int to_cmp (const void *restrict, const void *restrict);
static int size_cmp (const void *restrict, const void *restrict);

/* RMAT */

/* Recursively divide a grid of N x N by four to a single point, (i, j).
   Choose between the four quadrants with probability a, b, c, and d.
   Create an edge between node i and j.
*/

MTA ("mta expect parallel context") MTA ("mta inline")
static void
rmat_edge (int64_t * iout, int64_t * jout,
           int SCALE, double A, double B, double C, double D,
           const double *rn)
{
  static int xxx = 0;
  size_t rni = 0;
  int64_t i = 0, j = 0;
  int64_t bit = ((int64_t) 1) << (SCALE - 1);

  while (1) {
    const double r = rn[rni++];
    if (r > A) {                /* outside quadrant 1 */
      if (r <= A + B)           /* in quadrant 2 */
        j |= bit;
      else if (r <= A + B + C)  /* in quadrant 3 */
        i |= bit;
      else {                    /* in quadrant 4 */
        j |= bit;
        i |= bit;
      }
    }
    if (1 == bit)
      break;

    /*
      Assuming R is in (0, 1), 0.95 + 0.1 * R is in (0.95, 1.05).
      So the new probabilities are *not* the old +/- 10% but
      instead the old +/- 5%.
    */
    A *= (9.5 + rn[rni++]) / 10;
    B *= (9.5 + rn[rni++]) / 10;
    C *= (9.5 + rn[rni++]) / 10;
    D *= (9.5 + rn[rni++]) / 10;
    /* Used 5 random numbers. */

    {
      const double norm = 1.0 / (A + B + C + D);
      A *= norm;
      B *= norm;
      C *= norm;
    }
    /* So long as +/- are monotonic, ensure a+b+c+d <= 1.0 */
    D = 1.0 - (A + B + C);

    bit >>= 1;
  }
  /* Iterates SCALE times. */
  *iout = i;
  *jout = j;
}

MTA ("mta expect serial context")
void
rmat_edge_list (size_t listlen,
                int64_t * iout, size_t iout_stride,
                int64_t * jout, size_t jout_stride,
                int SCALE, double A, double B, double C, double D,
                size_t nwork, double *work,
                size_t niwork, int64_t * iwork)
{
  int64_t *restrict i = iout;
  int64_t *restrict j = jout;

  size_t nrand;
  double *restrict urand = work;
  size_t k, thislen;
  size_t *restrict which_self = (size_t *) iwork;
  int64_t *restrict tmpidx = &iwork[listlen + 1];

  /* Needs 5*SCALE random numbers per edge. */
  nrand = listlen * 5 * SCALE;
  assert (nrand <= nwork);
  assert (2 * listlen + 1 <= niwork);

  k = 0;
  while (k < listlen) {
    size_t thislen = listlen - k;
    size_t nselfedge = 0;
    int64_t *restrict ip = &i[k * iout_stride];
    int64_t *restrict jp = &j[k * iout_stride];

    assert (thislen * 5 * SCALE <= nwork);
    uniform_random (thislen * 5 * SCALE, urand);

    OMP ("omp parallel for")
      MTA ("mta assert nodep")
      MTASTREAMS ()
      for (size_t k2 = 0; k2 < thislen; ++k2) {
        rmat_edge (&ip[k2 * iout_stride], &jp[k2 * jout_stride],
                   SCALE, A, B, C, D, &urand[k2 * 5 * SCALE]);
        if (ip[k2 * iout_stride] == jp[k2 * jout_stride]) {
          size_t wi = stinger_size_fetch_add (&nselfedge, 1);
          which_self[wi] = k2;
        }
      }
    qsort (which_self, nselfedge, sizeof (*which_self), size_cmp);
    which_self[nselfedge] = thislen;

    if (nselfedge) {
      for (size_t k2 = 0; k2 < nselfedge; ++k2) {
        const size_t ntomove = which_self[k2 + 1] - which_self[k2];
        const size_t dest = which_self[k2];

        OMP ("omp parallel for")
          MTA ("mta assert nodep")
          MTASTREAMS ()
          for (size_t kidx = 0; kidx < ntomove; ++kidx)
            tmpidx[kidx] = ip[(dest + 1 + kidx) * iout_stride];
        OMP ("omp parallel for")
          MTA ("mta assert nodep")
          MTASTREAMS ()
          for (size_t kidx = 0; kidx < ntomove; ++kidx)
            ip[(dest + kidx) * iout_stride] = tmpidx[kidx];

        OMP ("omp parallel for")
          MTA ("mta assert nodep")
          MTASTREAMS ()
          for (size_t kidx = 0; kidx < ntomove; ++kidx)
            tmpidx[kidx] = jp[(dest + 1 + kidx) * jout_stride];
        OMP ("omp parallel for")
          MTA ("mta assert nodep")
          MTASTREAMS ()
          for (size_t kidx = 0; kidx < ntomove; ++kidx)
            jp[(dest + kidx) * jout_stride] = tmpidx[kidx];
      }
    }
    k = listlen - nselfedge;
  }
}

#define FROM_REF(k) (buf[2*(k)])
#define TO_REF(k) (buf[2*(k)+1])
#define FROM_REF2(k) (buf2[2*(k)])
#define TO_REF2(k) (buf2[2*(k)+1])
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

void
generate_undirected_rmat (int64_t ** off_out,
                          int64_t ** from_out, int64_t ** ind_out,
                          int64_t ** count_out, size_t n_undir_edge_gen,
                          int SCALE, double A, double B, double C, double D,
                          size_t nwork, double *work, size_t niwork,
                          int64_t * iwork, uint64_t min_size)
{
  const int64_t nv = (1u << SCALE);
  size_t ne, ne_gen;
  int64_t *restrict off = NULL;
  size_t *restrict off_generated = NULL;        /* edges may collapse; track offsets of expanded. */
  int64_t *restrict from = NULL;
  int64_t *restrict ind = NULL;
  int64_t *restrict count = NULL;

  int64_t *restrict buf_hndl = NULL;
  int64_t *restrict buf = NULL;
  int64_t *restrict buf2 = NULL;

  off = xcalloc (nv + 1, sizeof (*off));
  /* These next are overestimates, but not by much if the graph is sparse. */
  from = xmalloc (MAX(2 * n_undir_edge_gen, min_size) * sizeof (*ind));
  ind = xmalloc (MAX(2 * n_undir_edge_gen, min_size) * sizeof (*ind));
  count = xmalloc (MAX(2 * n_undir_edge_gen, min_size) * sizeof (*count));

  /* Remaining allocations are freed at the end. */
  /* Track offsets in the generated edges, not the compressed form. */
  off_generated = xcalloc (nv + 2, sizeof (*off));
  /* rows of <from, to>, doubled to hold reverse edges. */
  buf_hndl = xmalloc (2 * (4 * n_undir_edge_gen * sizeof (*buf_hndl)));
  buf = buf_hndl;
  buf2 = &buf[4 * n_undir_edge_gen];

  /* To fit, generate at most nwork / (5*scale) edges at a time. */
  ne_gen = nwork / (5 * SCALE);
  if (ne_gen < 1) {
    fprintf (stderr, "Too little workspace to generate even one edge.\n");
    abort ();
  }
  if (ne_gen > n_undir_edge_gen)
    ne_gen = n_undir_edge_gen;

  /* assert NOT parallel! */
  for (size_t ke = 0; ke < n_undir_edge_gen; ke += ne_gen) {
    const size_t len = (ke + ne_gen > n_undir_edge_gen ?
                        n_undir_edge_gen - ke : ne_gen);

    rmat_edge_list (len, &FROM_REF (2 * ke), 4, &TO_REF (2 * ke), 4,
                    SCALE, A, B, C, D, nwork, work, niwork, iwork);
  }

  /* Add reverse edges. */
  OMP ("omp parallel for")
    MTA ("mta assert nodep *buf")
    MTASTREAMS ()
    for (size_t ke = 0; ke < n_undir_edge_gen; ++ke) {
      FROM_REF (2 * ke + 1) = TO_REF (2 * ke);
      TO_REF (2 * ke + 1) = FROM_REF (2 * ke);
    }

  /* Initial histogram & cumulative sum for the outer sort */
  off_generated += 2;
  OMP ("omp parallel for")
    MTASTREAMS ()MTA ("mta assert nodep *buf")
    MTA ("mta assert nodep *off_generated")
    for (size_t ke = 0; ke < 2 * n_undir_edge_gen; ++ke)
      stinger_size_fetch_add (&off_generated[FROM_REF (ke)], 1);

  for (size_t i = 0; i < nv; ++i)
    off_generated[i] += off_generated[i - 1];
  off_generated--;

  /* Copy into place. */
  OMP ("omp parallel for")
    MTA ("mta assert nodep")
    MTASTREAMS ()
    for (size_t ke = 0; ke < 2 * n_undir_edge_gen; ++ke) {
      size_t dest = stinger_size_fetch_add (&off_generated[FROM_REF (ke)], 1);
      /* MTA won't parallelize: assert (dest < 2*n_undir_edge_gen); */
      FROM_REF2 (dest) = FROM_REF (ke);
      TO_REF2 (dest) = TO_REF (ke);
    }

  off_generated--;

  assert (off_generated[nv] = 2 * n_undir_edge_gen);
  /* Sort edge lists. */
  OMP ("omp parallel for")
    MTA ("mta assert parallel")
    MTASTREAMS ()MTA ("mta dynamic schedule")
    for (size_t k = 0; k < nv; ++k) {
      const size_t deg = off_generated[k + 1] - off_generated[k];
      qsort (&FROM_REF2 (off_generated[k]), deg, 2 * sizeof (*buf2), to_cmp);
    }

  buf = buf2;

  /* Histogram unique edges. */
  if (FROM_REF (0) != TO_REF (0))
    off[FROM_REF (0)] = 1;
  else
    off[FROM_REF (0)] = 0;

  OMP ("omp parallel for")
    MTASTREAMS ()MTA ("mta assert nodep *buf")
    MTA ("mta assert nodep *off")
    for (size_t ke = 1; ke < 2 * n_undir_edge_gen; ++ke) {
      if (!(FROM_REF (ke) == FROM_REF (ke - 1)
            && TO_REF (ke) == TO_REF (ke - 1)) && FROM_REF (ke) != TO_REF (ke)) {
#if !defined(__MTA__)
        /* Yes, their compiler is this stupid. */
        stinger_int64_fetch_add (&off[FROM_REF (ke)], 1);
#else
        int_fetch_add (&off[FROM_REF (ke)], 1);
#endif
      }
    }

  /* Cumulative sum & placement */
  ne = 0;
  MTASTREAMS ()
    for (size_t i = 0; i < nv; ++i) {
      size_t tmp = off[i];
      off[i] = ne;
      ne += tmp;
    }
  off[nv] = ne;

  /* XXX: Could realloc ind, count (in that order) here. */

  /* Copy & count. */
  OMP ("omp parallel for")
    MTA ("mta assert nodep *off")
    MTA ("mta assert nodep *off_generated")
    MTA ("mta assert nodep *from")
    MTA ("mta assert nodep *ind")
    MTA ("mta assert nodep *count")
    MTASTREAMS ()MTA ("mta dynamic schedule")
    for (size_t i = 0; i < nv; ++i) {
      int64_t prev_j = -1;
      size_t dst = off[i];
      const size_t src_end = off_generated[i + 1];
      /* assert NOT parallel to maintain the sorted order. */
      for (size_t src = off_generated[i]; src < src_end; ++src) {
        const int64_t j = TO_REF (src);
        if (j != i) {
          if (j != prev_j) {
            from[dst] = i;
            ind[dst] = j;
            count[dst] = 1;
            ++dst;
            prev_j = j;
          } else
            ++count[dst - 1];
        }
      }
    }

  *from_out = from;
  *off_out = off;
  *ind_out = ind;
  *count_out = count;

  /* Free workspace. */
  free (buf_hndl);
  free (off_generated);
}

int
to_cmp (const void *restrict avp, const void *restrict bvp)
{
  const int64_t *restrict ap = (int64_t *) avp;
  const int64_t *restrict bp = (int64_t *) bvp;
  return ap[1] - bp[1];
}

int
size_cmp (const void *restrict avp, const void *restrict bvp)
{
  const size_t a = *(size_t *) avp;
  const size_t b = *(size_t *) bvp;
  if (a > b)
    return 1;
  if (a < b)
    return -1;
  return 0;
}

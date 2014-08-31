/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(__MTA__)
#include <machine/runtime.h>
#include <mta_rng.h>
#endif

#include "stinger.h"
#include "xmalloc.h"


/* PRNG helpers */
#if !defined(__MTA__)
#if !defined(NRNG)
#define NRNG 8
#endif
static int nmtrng;
static unsigned short *eseed;
#else /* __MTA__ */
#if !defined(NRNG)
#define NRNG 100
#endif
static sync int only_one$ = 1;
static int nmtrng;
static RNG *restrict mtrng;
#endif

void
init_urandom (void)
{
  static long seed = 39362;
  const char *seedstr = NULL;

  seedstr = getenv ("SEED");
  if (seedstr)
    seed = strtol (seedstr, NULL, 10);
#if !defined(__MTA__)
  srand48 (seed);
#if defined(_OPENMP)
  nmtrng = omp_get_max_threads ();
#else
  nmtrng = 1;
#endif
  if (nmtrng < NRNG)
    nmtrng = NRNG;
  eseed = xmalloc (3 * nmtrng * sizeof (*eseed));
  for (size_t k = 0; k < NRNG * 3; ++k)
    eseed[k] = (unsigned short) lrand48 ();
  OMP ("omp barrier")
#else
    nmtrng = MTA_NUM_STREAMS ();  /* which can change before we generate. sigh. */
  if (nmtrng < NRNG)
    nmtrng = NRNG;
  mtrng = init_rn (nmtrng, seed);
#endif
}

void
free_urandom (void)
{
#if !defined(__MTA__)
  free (eseed);
#endif
}

MTA ("mta expect serial context")
void uniform_random (size_t nrand, double *rn)
{
#if !defined(__MTA__)
  OMP ("omp parallel for")
    for (int T = 0; T < nmtrng; ++T) {
      const size_t len = nrand / nmtrng;
      const size_t start = len * T;
      const size_t finish = len * (T + 1);
      for (size_t k = start; k < finish; ++k)
        do {
          assert (k < nrand);
          rn[k] = erand48 (&eseed[3 * T]);
        }
        while (rn[k] == 0.0 || rn[k] == 1.0);     /* Exclude endpoints. */
    }
#else
  int x = only_one$;
  MTA ("mta assert nodep")
    MTA ("mta assert noalias *rn")
    MTASTREAMS ()
    for (int T = 0; T < nmtrng; ++T) {
      const size_t len = nrand / nmtrng;
      const size_t start = len * T;
      const size_t finish = len * (T + 1);
      for (size_t k = start; k < finish; ++k)
        do {
          rn[k] = get_random_dbl (&mtrng[T]);     /* XXX: Docs say in [0,10), experiments say (0,1]. */
        }
        while (rn[k] == 0.0 || rn[k] == 1.0);     /* Exclude endpoints. */
    }
  only_one$ = x;
#endif
}

MTA ("mta expect parallel context")
double
poisson (const double lambda, const double *rn, size_t * rni,
         const size_t nrand)
/* Knuth's algorithm to generate a Poisson-distributed random number
   by simulation.  This is not the fastest method, but it works
   directly from uniform variates in (0, 1).  Will need O(lambda)
   random numbers.  */
{
  const double L = exp (-lambda);
  size_t i = *rni;
  int k = 0;
  double p = 1.0;
  do {
    if (i >= nrand)
      return -1;                /* Ran out of random numbers. */
    p *= rn[i++];
    ++k;
  }
  while (p > L);
  *rni = i;
  return k - 1;
}

MTA ("mta expect parallel context")
double exponential (const double lambda, const double r)
{
  return -log (r) / lambda;
}


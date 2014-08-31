/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#include "stinger-atomics.h"
#include "stinger-utils.h"
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

void count_all_triangles (const size_t nv, const int64_t * off, const int64_t * ind, int64_t * ntri);

static int
i64cmp (const void *ap, const void *bp)
{
  return (int) ((*(int64_t *) ap) - (*(int64_t *) bp));
}
int csr_cmp(const void * a, const void * b);

#define VERIFY 0

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

  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, -2);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush(stdout);

  int64_t * ntri = xcalloc (nv, sizeof (*ntri));
  int64_t * local_cc = xcalloc (nv, sizeof (*local_cc));
  int64_t * affected = xcalloc (nv, sizeof (*affected));
  int64_t * result = xcalloc (batch_size, sizeof (*result));
  double * update_time_trace = xcalloc (nbatch, sizeof(*update_time_trace));

  /* Calculate the initial clustering coefficients */
  tic();
  count_all_triangles (nv, off, ind, ntri);
  PRINT_STAT_DOUBLE ("time_triangles", toc());
  fflush(stdout);

  tic();
  int64_t global_degsum = 0;
  int64_t global_ntri = 0;
  OMP("omp parallel for")
  MTA("mta assert nodep")
  for (int64_t i = 0; i < nv; i++) {
    int64_t deg = off[i+1] - off[i];
    int64_t d = deg * (deg-1);
    local_cc[i] = (d ? ntri[i] / (double) d : 0.0);
    global_degsum += d;
    global_ntri += ntri[i];
  }
  double global_cc = (global_degsum ? global_ntri / (double) global_degsum : 0.0);
  PRINT_STAT_INT64 ("global_ntri_init", global_ntri);
  PRINT_STAT_INT64 ("global_degsum_init", global_degsum);
  PRINT_STAT_DOUBLE ("initial_cc", global_cc);
  PRINT_STAT_DOUBLE ("time_initial_cc", toc());
  PRINT_STAT_HEX64("error_code", stinger_consistency_check (S, nv));
  fflush(stdout);

  /* Updates */
  int64_t ntrace = 0;

  for (int64_t actno = 0; actno < nbatch * batch_size; actno += batch_size)
  {
    tic();

    const int64_t endact = (actno + batch_size > naction ? naction : actno + batch_size);

    int64_t *actions = &action[2*actno];
    int64_t numActions = endact - actno;

    OMP("omp parallel for")
    MTA("mta assert parallel")
    for (int64_t k = 0; k < nv; k++)
    {
      affected[k] = 0;
    }

    OMP("omp parallel for")
    MTA("mta assert parallel")
    for (int64_t k = 0; k < batch_size; k++)
    {
      result[k] = 0;
    }

    /* For each unique vertex represented in the batch, we remove its degree contribution
       from the global degree sum, or the denominator of the global clustering coefficient. */

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    OMP("omp parallel for")
    for (uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];

      if (i != j && i >= 0) {
	int64_t r1 = stinger_insert_edge (S, 0, i, j, 1, actno+1);
	int64_t r2 = stinger_insert_edge (S, 0, j, i, 1, actno+1);
	stinger_int64_fetch_add (&affected[i], 1);
	stinger_int64_fetch_add (&affected[j], 1);

	if (r1 == 1 && r2 == 1) result[k] = 3;
	else if (r1 == 1) result[k] = 1;
	else if (r2 == 1) result[k] = 2;
      }
    }

    MTA("mta assert parallel")
    MTA("mta block dynamic schedule")
    for (uint64_t k = 0; k < endact - actno; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];

      if (result[k] == 1 && i > j) continue;

      if (i != j && i >= 0 && result[k]) {
	int64_t i_deg, j_deg;
	i_deg = stinger_outdegree(S, i);

	int64_t * i_neighbors = xmalloc (i_deg * sizeof(*i_neighbors));
	int64_t * i_neighbors_ts = xmalloc (i_deg * sizeof(*i_neighbors_ts));
	size_t i_degree;
	stinger_gather_successors (S, i, &i_degree, i_neighbors, NULL, i_neighbors_ts, NULL, NULL, i_deg);
	assert (i_degree == i_deg);

	int64_t ** pointers = xmalloc(i_degree * sizeof(int64_t *));
	for(int64_t j = 0; j < i_degree; j++) {
	  pointers[j] = i_neighbors+j; 
	}

	qsort (pointers, i_degree, sizeof(int64_t *), csr_cmp);

	int64_t * tmp_buffer = xmalloc(i_degree * sizeof(int64_t));
	for(int64_t j = 0; j < i_degree; j++) {
	  tmp_buffer[j] = *pointers[j];
	}
	for(int64_t j = 0; j < i_degree; j++) {
	  i_neighbors[j] = tmp_buffer[j];
	}
	for(int64_t i = 0; i < i_degree; i++) {
	  tmp_buffer[i] = i_neighbors_ts[pointers[i] - i_neighbors];
	}
	for(int64_t i = 0; i < i_degree; i++) {
	  i_neighbors_ts[i] = tmp_buffer[i];
	}

	int64_t incr = 0;

	j_deg = stinger_outdegree(S, j);
	size_t j_degree;

	if (!j_deg)
	  continue;

	int64_t * j_neighbors = xmalloc (j_deg * sizeof(*j_neighbors));
	int64_t * j_neighbors_ts = xmalloc (j_deg * sizeof(*j_neighbors_ts));
	stinger_gather_successors (S, j, &j_degree, j_neighbors, NULL, j_neighbors_ts, NULL, NULL, j_deg);
	assert (j_degree == j_deg);

	for (int64_t k2 = 0; k2 < j_degree; k2++)
	{
	  int64_t u = j_neighbors[k2];
	  int64_t where = find_in_sorted (u, i_degree, i_neighbors);
	  if (where >= 0) {
	    int64_t neighbor_i_ts = i_neighbors_ts[where];
	    int64_t neighbor_j_ts = j_neighbors_ts[k2];
	    if (neighbor_i_ts != actno+1 && neighbor_j_ts != actno+1) {
	      incr += 1;
	      stinger_int64_fetch_add (&ntri[u], 2);
	      stinger_int64_fetch_add (&affected[u], 1);
	    } else if (i < u && neighbor_i_ts != actno+1) {
	      incr += 1;
	      stinger_int64_fetch_add (&ntri[u], 2);
	      stinger_int64_fetch_add (&affected[u], 1);
	    } else if (j < u && neighbor_j_ts != actno+1) {
	      incr += 1;
	      stinger_int64_fetch_add (&ntri[u], 2);
	      stinger_int64_fetch_add (&affected[u], 1);
	    } else if (neighbor_i_ts == actno+1 && neighbor_j_ts == actno+1) {
	      if (i < u && j < u) {
		incr += 1;
		stinger_int64_fetch_add (&ntri[u], 2);
		stinger_int64_fetch_add (&affected[u], 1);
	      }
	    }
	  }
	}

	free (j_neighbors_ts);
	free (j_neighbors);
	free (pointers);
	free (i_neighbors_ts);
	free (i_neighbors);

	stinger_int64_fetch_add (&ntri[i], 2*incr);
	stinger_int64_fetch_add (&ntri[j], 2*incr);
	global_ntri += (6*incr);
      }
    }

    global_degsum = 0;
    OMP("omp parallel for")
    MTA("mta assert nodep")
    for (int64_t i = 0; i < nv; i++) {
      int64_t deg = stinger_outdegree (S, i);
      int64_t d = deg * (deg-1);
      global_degsum += d;
    }

    /* Update the local clustering coefficient for affected vertices */
    OMP("omp parallel for")
    MTA("mta assert parallel")
    for (int64_t k = 0; k < nv; k++) {
      if (affected[k]) {
	double d = (double) stinger_outdegree (S, k);
	double new_local_cc = (d > 1 ? ntri[k] / (d * (d-1)) : 0.0);
	local_cc[k] = new_local_cc;
      }
    }

    update_time_trace[ntrace] = toc();
    ++ntrace;
  } /* End of batch */

  /* Print the times */
  double time_updates = 0;
  for (int64_t k = 0; k < nbatch; k++) {
    time_updates += update_time_trace[k];
  }
  global_cc = (global_degsum ? global_ntri / (double) global_degsum : 0.0);
  PRINT_STAT_DOUBLE ("final_cc", global_cc);
  PRINT_STAT_DOUBLE ("time_updates", time_updates);
  PRINT_STAT_DOUBLE ("updates_per_sec", (nbatch * batch_size) / time_updates); 

  PRINT_STAT_INT64 ("global_ntri_final", global_ntri);
  PRINT_STAT_INT64 ("global_degsum_final", global_degsum);

  PRINT_STAT_HEX64("error_code", stinger_consistency_check (S, nv));

#if VERIFY
  /* Calculate the correct clustering coefficients */
  tic();
  int64_t * ntri_final = xcalloc (nv, sizeof (*ntri));
  stinger_to_sorted_csr (S, nv, &off, &ind, NULL, NULL, NULL, NULL);
  count_all_triangles (nv, off, ind, ntri_final);
  global_degsum = 0;
  global_ntri = 0;
  OMP("omp parallel for")
  MTA("mta assert nodep")
  for (int64_t i = 0; i < nv; i++) {
    int64_t deg = off[i+1] - off[i];
    int64_t d = deg * (deg-1);
    local_cc[i] = (d ? ntri_final[i] / (double) d : 0.0);
    global_degsum += d;
    global_ntri += ntri_final[i];
  }
  PRINT_STAT_DOUBLE ("time_verify", toc());
  PRINT_STAT_INT64 ("global_ntri_verify", global_ntri);
  PRINT_STAT_INT64 ("global_degsum_verify", global_degsum);

  int64_t total_errors = 0;
  MTA("mta assert nodep")
  for (int64_t i = 0; i < nv; i++) {
    if (ntri_final[i] != ntri[i]) {
      total_errors += 1;
      printf("\n(%ld),", i);
    }
  }
  PRINT_STAT_INT64 ("total_errors", total_errors);
#endif

  fflush(stdout);

  free (update_time_trace);
  free (affected);
  free (local_cc);
  free (ntri);

  stinger_free_all (S);
  free (graphmem);
  free (actionmem);
  STATS_END();
}



MTA("mta expect parallel context") MTA("mta inline")
size_t
count_intersections (const int64_t ai,
		     const size_t alen,
		     const int64_t * a,
		     const int64_t bi,
		     const size_t blen,
		     const int64_t * b)
{
  size_t ka = 0, kb = 0;
  size_t out = 0;

  if (!alen || !blen || a[alen-1] < b[0] || b[blen-1] < a[0])
    return 0;

  while (1) {
    if (ka >= alen || kb >= blen) break;

    int64_t va = a[ka];
    int64_t vb = b[kb];

    /* Skip self-edges. */
    if (UNLIKELY(va == ai)) {
      ++ka;
      if (ka >= alen) break;
      va = a[ka];
    }
    if (UNLIKELY(vb == bi)) {
      ++kb;
      if (kb >= blen) break;
      vb = b[kb];
    }

    if (va == vb) { ++ka; ++kb; ++out; }
    else if (va < vb) {
      /* Advance ka */
      ++ka;
      while (ka < alen && a[ka] < vb) ++ka;
    } else {
      /* Advance kb */
      ++kb;
      while (kb < blen && va > b[kb]) ++kb;
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
  {
    const size_t N = nv;
    OMP ("omp for schedule(dynamic,128)")
    MTA ("mta dynamic schedule") MTASTREAMS ()
    for (size_t i = 0; i < N; ++i)
      ntri[i] = count_triangles (nv, off, ind, i);
  }
}

/* -*- mode: C; mode: folding; fill-column: 70; -*- */

#include "stinger.h"
#include "prng.h"
#include "timer.h"
#include "xmalloc.h"

/* In rmat.c */
extern void rmat_edge_list (size_t listlen,
                            int64_t * iout, size_t iout_stride,
                            int64_t * jout, size_t jout_stride,
                            int SCALE, double A, double B, double C, double D,
                            size_t nwork, double *work,
                            size_t niwork, int64_t * iwork);
extern void generate_undirected_rmat (int64_t ** off_out,
                                      int64_t ** from_out, int64_t ** ind_out,
                                      int64_t ** count_out,
                                      size_t n_undir_edge_gen,
                                      int SCALE, double A, double B,
                                      double C, double D,
                                      size_t nwork, double *work,
                                      size_t niwork, int64_t * iwork, uint64_t min_size);


/* {{{ Parameters and defaults */

/* Number of additions/deletions to run. */
#define NEDGE_ACTIONS_DEFAULT 2000
static int nedge_actions = NEDGE_ACTIONS_DEFAULT;

/* Scale parameter for RMAT(scale), 2**scale vertices */
//#define SCALE_DEFAULT 20
#define SCALE_DEFAULT 16
#define SCALE_MAX 40
static int scale = SCALE_DEFAULT;

/* Edge multiplier, edgefact * 2**scale edges in generated graph. */
#define EDGEFACT_DEFAULT 8
static int edgefact = EDGEFACT_DEFAULT;

/* RMAT probabilities for the quadrants [A, B; C, D]. */
#define PROB_A_DEFAULT 0.55
#define PROB_B_DEFAULT 0.1
#define PROB_C_DEFAULT 0.1
#define PROB_D_DEFAULT 0.25
static double prob_a = PROB_A_DEFAULT, prob_b = PROB_B_DEFAULT,
  prob_c = PROB_C_DEFAULT, prob_d = PROB_D_DEFAULT;

/* Probability that an action will choose from the deletion queue (from U(0,1)). */
#define PROB_DELETE_DEFAULT 0.0625
static double prob_delete = PROB_DELETE_DEFAULT;
/* Probability of appending an edge to the deletion queue (from U(0,1)). */
#define PROB_DELETEQ_DEFAULT 0.0625
static double prob_deleteq = PROB_DELETEQ_DEFAULT;

#define INITIAL_GRAPH_NAME_DEFAULT "initial-graph.bin"
static const char *initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
#define ACTION_STREAM_NAME_DEFAULT "action-stream.bin"
static const char *action_stream_name = ACTION_STREAM_NAME_DEFAULT;

/* }}} */

/* {{{ Command line parsing */

static void
usage (FILE * out, char *progname)
{
  fprintf (out,
           "%s [--scale=SCALE] [--stats=A,B,C,D] [--edgefact=EDGEFACT] [--del=DELETE] [--delq=DELETEQ]\n"
           "\t\t[--nact=NACT] [initial-graph.bin [action-stream.bin]]\n"
           "\tDefaults:\n"
           "\t   NACT = %d\n"
           "\t   SCALE = %d, EDGEFACT = %d,\n"
           "\t   A = %g, B = %g, C = %g, D = %g\n"
           "\t   DELETE = %g, DELETEQ = %g\n"
           "\t   initial-graph name = \"%s\"\n"
           "\t   action-stream name = \"%s\"\n",
#if !defined(__MTA__)
           basename (progname),
#else
           progname,
#endif
           NEDGE_ACTIONS_DEFAULT,
           SCALE_DEFAULT, EDGEFACT_DEFAULT,
           PROB_A_DEFAULT, PROB_B_DEFAULT, PROB_C_DEFAULT, PROB_D_DEFAULT,
           (double) PROB_DELETE_DEFAULT, (double) PROB_DELETEQ_DEFAULT,
           INITIAL_GRAPH_NAME_DEFAULT, ACTION_STREAM_NAME_DEFAULT);
}

static void
parse_args (const int argc, char *argv[])
{
  int arg_err = 0;
  while (1) {
    int c;
    static struct option long_options[] = {
      {"scale", required_argument, 0, 's'},
      {"edgefact", required_argument, 0, 'E'},
      {"nact", required_argument, 0, 'a'},
      {"del", required_argument, 0, 'd'},
      {"delq", required_argument, 0, 'D'},
      {"stats", required_argument, 0, 'S'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    extern char *optarg;
    extern int optind, opterr, optopt;

    c = getopt_long (argc, argv, "s:E:a:d:D:S:h?",
                     long_options, &option_index);
    if (-1 == c)
      break;
    switch (c) {
    case '?':
    case 'h':
      usage (stdout, argv[0]);
      exit (EXIT_SUCCESS);
    case 's':
      errno = 0;
      scale = strtol (optarg, NULL, 10);
      if (errno || scale < 0 || scale > SCALE_MAX) {
        fprintf (stderr,
                 "Error: scale outside range {0, 1, ..., %d}: %s\n",
                 SCALE_MAX, optarg);
        arg_err = 1;
      }
      break;
    case 'E':
      errno = 0;
      edgefact = strtol (optarg, NULL, 10);
      if (errno || edgefact < 1) {
        fprintf (stderr,
                 "Error: edgefact outside range {1, 2, ..., ∞}: %s\n",
                 optarg);
        arg_err = 1;
      }
      break;
    case 'a':
      errno = 0;
      nedge_actions = strtol (optarg, NULL, 10);
      if (errno || nedge_actions < 0) {
        fprintf (stderr,
                 "Error: number of actions outside range {0, 1, ..., ∞}: %s\n",
                 optarg);
        arg_err = 1;
      }
      break;
    case 'd':
      errno = 0;
      prob_delete = strtod (optarg, NULL);
      if (errno || prob_delete < 0 || prob_delete > 1.0) {
        fprintf (stderr,
                 "Error: delete probability outside range [0, 1]: %s\n",
                 optarg);
        arg_err = 1;
      }
      break;
    case 'D':
      errno = 0;
      prob_deleteq = strtod (optarg, NULL);
      if (errno || prob_deleteq < 0 || prob_deleteq > 1.0) {
        fprintf (stderr,
                 "Error: delete enqueue probability outside range [0, 1]: %s\n",
                 optarg);
        arg_err = 1;
      }
      break;
    case 'S':
      {
        char *subarg = NULL;
        char *subarg_next = NULL;
        errno = 0;
        subarg = optarg;
        prob_a = strtod (subarg, &subarg_next);
        if (subarg == subarg_next)
          prob_a = -1.0;
        subarg = subarg_next;
        if (!isspace (*subarg))
          ++subarg;
        prob_b = strtod (subarg, &subarg_next);
        if (subarg == subarg_next)
          prob_b = -1.0;
        subarg = subarg_next;
        if (!isspace (*subarg))
          ++subarg;
        prob_c = strtod (subarg, &subarg_next);
        if (subarg == subarg_next)
          prob_c = -1.0;
        subarg = subarg_next;
        if (!isspace (*subarg))
          ++subarg;
        prob_d = strtod (subarg, &subarg_next);
        if (subarg == subarg_next)
          prob_d = -1.0;
        subarg = subarg_next;
        if (errno) {
          fprintf (stderr, "Error parsing probability string %s\n", optarg);
          arg_err = 1;
        } else {
          long double sum = 0;
          int nempty = 0;
          if (prob_a == -1)
            ++nempty;
          else if (prob_a < 0 || prob_a > 1) {
            fprintf (stderr, "Error: A outside range [0, 1]: %s\n", optarg);
            arg_err = 1;
          } else
            sum += prob_a;
          if (prob_b == -1)
            ++nempty;
          else if (prob_b < 0 || prob_b > 1) {
            fprintf (stderr, "Error: B outside range [0, 1]: %s\n", optarg);
            arg_err = 1;
          } else
            sum += prob_b;
          if (prob_c == -1)
            ++nempty;
          else if (prob_c < 0 || prob_c > 1) {
            fprintf (stderr, "Error: C outside range [0, 1]: %s\n", optarg);
            arg_err = 1;
          } else
            sum += prob_c;
          if (prob_d == -1)
            ++nempty;
          else if (prob_d < 0 || prob_d > 1) {
            fprintf (stderr, "Error: D outside range [0, 1]: %s\n", optarg);
            arg_err = 1;
          } else
            sum += prob_d;
          sum = (1.0 - sum) / nempty;
          if (prob_a == -1)
            prob_a = sum;
          if (prob_b == -1)
            prob_b = sum;
          if (prob_c == -1)
            prob_c = sum;
          if (prob_d == -1)
            prob_d = sum;
        }
        break;
      }
    default:
      arg_err = 1;
      break;
    }
  }
  if (arg_err) {
    usage (stderr, argv[0]);
    exit (EXIT_FAILURE);
  }

  if (optind < argc)
    initial_graph_name = argv[optind++];
  if (optind < argc)
    action_stream_name = argv[optind++];
}

/* }}} */

/* {{{ Output files */

#if !defined(__MTA__)
static int rmat_fd, action_fd;
#endif

static void
open_output_files (void)
{
#if !defined(__MTA__)
  mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

  rmat_fd = open (initial_graph_name, O_WRONLY | O_CREAT | O_TRUNC, mode);
  if (rmat_fd < 0) {
    perror ("Couldn't open output graph: ");
    abort ();
  }

  action_fd = open (action_stream_name, O_WRONLY | O_CREAT | O_TRUNC, mode);
  if (action_fd < 0) {
    perror ("Couldn't open action stream: ");
    abort ();
  }
#else /* __MTA__ inanity */
  extern void xmt_luc_io_init (void);
  xmt_luc_io_init ();
#endif
}

#if !defined(__MTA__)
static ssize_t
xwrite (int fd, const void *buf, size_t len)
{
  ssize_t lenwritten;
  ssize_t total_bytes = 0;
  while (len) {
    lenwritten = write (fd, buf, len);
    if (lenwritten < 0)
      return lenwritten;
    buf += lenwritten;
    len -= lenwritten;
    total_bytes += lenwritten;
  }
  return total_bytes;
}
#endif

static void
dump_initial_graph (int64_t nv, int64_t ne, const int64_t * restrict off,
                    const int64_t * restrict from,
                    const int64_t * restrict ind,
                    const int64_t * restrict weight, const size_t niwork,
                    int64_t * restrict iwork)
{
  const uint64_t endian_check = 0x1234ABCDul;

#if !defined(__MTA__)
  /* The sane option. */
  if (sizeof (endian_check) != xwrite (rmat_fd, &endian_check,
                                       sizeof (endian_check)))
    goto err;

  if (sizeof (nv) != xwrite (rmat_fd, &nv, sizeof (nv)))
    goto err;

  if (sizeof (ne) != xwrite (rmat_fd, &ne, sizeof (ne)))
    goto err;

  if ((nv + 1) * sizeof (*off) !=
      xwrite (rmat_fd, off, (nv + 1) * sizeof (*off)))
    goto err;

  /* Unused.
     if (ne * sizeof (*from) != xwrite (rmat_fd, from, ne * sizeof (*from)))
     goto err;
  */

  if (ne * sizeof (*ind) != xwrite (rmat_fd, ind, ne * sizeof (*ind)))
    goto err;

  if (ne * sizeof (*weight) !=
      xwrite (rmat_fd, weight, ne * sizeof (*weight)))
    goto err;
#else /* __MTA__ inanity */
  extern void xmt_luc_snapout (const char *, void *, size_t);
  assert (niwork >= 3 + (nv + 1) + 2 * ne);
  iwork[0] = endian_check;
  iwork[1] = nv;
  iwork[2] = ne;
  for (size_t k = 0; k <= nv; ++k)
    iwork[2 + k + 1] = off[k];
  for (size_t k = 0; k < ne; ++k)
    iwork[2 + (nv + 1) + k + 1] = ind[k];
  for (size_t k = 0; k < ne; ++k)
    iwork[2 + (nv + 1) + ne + k + 1] = weight[k];
  xmt_luc_snapout (initial_graph_name, iwork,
                   (3 + (nv + 1) + 2 * ne) * sizeof (*iwork));
#endif

  return;
 err:
  perror ("Error writing initial rmat graph: ");
  abort ();
}

static void
dump_actions (int64_t naction, const int64_t * restrict action,
              const size_t niwork, int64_t * restrict iwork)
{
  const uint64_t endian_check = 0x1234ABCDul;

#if !defined(__MTA__)
  if (sizeof (endian_check) != xwrite (action_fd, &endian_check,
                                       sizeof (endian_check)))
    goto err;

  if (sizeof (naction) != xwrite (action_fd, &naction, sizeof (naction)))
    goto err;

  if (2 * naction * sizeof (*action) != xwrite (action_fd, action,
                                                2 * naction *
                                                sizeof (*action)))
    goto err;
#else /* ___MTA__ */
  extern void xmt_luc_snapout (const char *, void *, size_t);
  assert (niwork >= 2 + 2 * naction);
  iwork[0] = endian_check;
  iwork[1] = naction;
  for (size_t k = 0; k < 2 * naction; ++k)
    iwork[1 + k + 1] = action[k];
  xmt_luc_snapout (action_stream_name, iwork,
                   (2 + 2 * naction) * sizeof (*iwork));
#endif

  return;
 err:
  perror ("Error writing action list: ");
  abort ();
}

static void
close_output (void)
{
#if !defined(__MTA__)
  if (action_fd >= 0)
    close (action_fd);
  if (rmat_fd >= 0)
    close (rmat_fd);
#endif
}

/* }}} */

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

static int64_t nv, ne, *off, *from, *ind, *weight;
static int64_t *elist;
static int64_t ndelq, *delq;
static int edge_k, del_k;
static int64_t *actions;

static size_t nrnd;
static double *rnd;

static size_t niwork;
static int64_t *iwork;

static double rmat_time;

int
main (int argc, char **argv)
{
  parse_args (argc, argv);
  open_output_files ();
  init_urandom ();
  init_timer ();

  nv = ((int64_t) 1) << scale;
  printf ("nv: %d\n", (int) nv);

  /* scratch space for random number generation */
  nrnd = 5 * scale * (2 * nv);
  if (nrnd < 5 * scale * nedge_actions)
    nrnd = 5 * scale * nedge_actions;
  printf ("nrnd: %d\n", (int) nrnd);
  fflush(stdout);
  rnd = xmalloc (nrnd * sizeof (*rnd));

  /* scratch space for removing duplicates, packing output (on the bone-headed machine) */
#if !defined(__MTA__)
  niwork = MAX(2 * edgefact * nv + 1, 2 * nedge_actions+ 1);
#else /* __MTA__ */
  niwork = 3 + (nv + 1) + 4 * edgefact * nv;
#endif
  iwork = xmalloc (niwork * sizeof (*iwork));

  /* Generate static RMAT */
  tic ();
  generate_undirected_rmat (&off, &from, &ind, &weight, nv * edgefact,
                            scale, prob_a, prob_b, prob_c, prob_d,
                            nrnd, rnd, niwork, iwork, 2 * nedge_actions);
  ne = off[nv];
  rmat_time = toc ();
  printf ("time rmat: %g\n", rmat_time);
  printf ("ne: %d\n", (int) ne);
  fflush(stdout);

  dump_initial_graph (nv, ne, off, from, ind, weight, niwork, iwork);

  /* Save some edges as deletions */
  elist = xmalloc (ne * sizeof (*elist));
  uniform_random (nedge_actions, rnd);

  MTA ("mta assert parallel")
    for (int64_t k = 0; k < ne; k++) {
      elist[k] = k;
    }

  MTA ("mta assert parallel")
    for (int64_t k = 0; k < MIN(nedge_actions, ne); ++k) {
      const int64_t nleft = ne - k;
#if defined(__MTA__)
#define lrint(x) rint((x))
#endif
      const int64_t which = k + lrint (nleft * rnd[k]);
      const int64_t tmp = elist[which];
      elist[which] = elist[k];
      elist[k] = tmp;
    }
  uniform_random (5 * scale * nedge_actions, rnd);
  delq = xmalloc (2 * nedge_actions * sizeof (*delq));
  ndelq = 0;
  for (int64_t k = 0;
       k < ne && k < 5 * scale * nedge_actions && ndelq < nedge_actions; ++k)
    {
      if (rnd[k] < prob_deleteq) {
        delq[2 * ndelq] = from[elist[k]];
        delq[2 * ndelq + 1] = ind[elist[k]];
        ++ndelq;
      }
    }

  /* Generate the action list */
  nrnd = 5 * scale * nedge_actions;
  uniform_random (nrnd, rnd);
  rmat_edge_list (nedge_actions, from, 1, ind, 1,
                  scale, prob_a, prob_b, prob_c, prob_d, nrnd, rnd,
                  niwork, iwork);
  actions = xmalloc (2 * nedge_actions * sizeof (*actions));
  uniform_random (nedge_actions, rnd);
  edge_k = 0;
  del_k = 0;
  for (int k = 0; k < nedge_actions; ++k) {
    double flip = rnd[k];
    if (flip >= prob_delete) {  /* Insert a new edge */
      actions[2 * k] = from[edge_k];
      actions[2 * k + 1] = ind[edge_k];
      ++edge_k;
    } else {                      /* Delete an edge. */
      /* This is a bit nasty.  Choose from the queue of *old* edges
         until we run out, then pick from generated edges. */
      if (ndelq) {
        /* Use the queue from initial edges.  It was permuted, so
           pop off the end. */
        --ndelq;
        actions[2 * k] = -(delq[2 * ndelq] + 1);
        actions[2 * k + 1] = -(delq[2 * ndelq + 1] + 1);
      } else if (del_k < k) {
        /* Pick from added edges. */
        actions[2 * k] = -(actions[2 * del_k] + 1);
        actions[2 * k + 1] = -(actions[2 * del_k + 1] + 1);
        ++del_k;
      } else {
        /* This *could* happen, but the chances are so remote that
           it's not worth coping. */
        fprintf (stderr,
                 "The asoundingly unlikely happened: deletions "
                 "caught up to insertions.\n");
        abort ();
      }
    }
  }

  dump_actions (nedge_actions, actions, niwork, iwork);

  close_output ();

  free_urandom ();

  free (actions);
  free (delq);
  free (elist);
  free (weight);
  free (ind);
  free (from);
  free (off);
  free (iwork);
  free (rnd);

  return EXIT_SUCCESS;
}

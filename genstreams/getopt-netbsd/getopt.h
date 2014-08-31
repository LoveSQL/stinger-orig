#if !defined(GETOPT_HEADER_)
#define GETOPT_HEADER_

#define no_argument        0
#define required_argument  1
#define optional_argument  2

int getopt(int argc, char * const *, const char *);
extern char *optarg;
extern int optind, opterr, optopt;

struct option {
  const char *name;
  int has_arg;
  int *flag;
  int val;
};

int getopt_long(int, char * const *, const char *, const struct option *,
		int *);

#endif /* GETOPT_HEADER_ */

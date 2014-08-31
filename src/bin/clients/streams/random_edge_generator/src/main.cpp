#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/types.h>
#include <time.h>
#include <netdb.h>

extern "C" {
#include "stinger_core/stinger.h"
#include "stinger_core/xmalloc.h"
#include "stinger_utils/csv.h"
#include "stinger_utils/timer.h"
#include "stinger_utils/stinger_sockets.h"
}

#include "random_edge_generator.h"
#include "build_name.h"


using namespace gt::stinger;

#define E_A(X,...) fprintf(stderr, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define E(X) E_A(X,NULL)
#define V_A(X,...) fprintf(stdout, "%s %s %d:\n\t" #X "\n", __FILE__, __func__, __LINE__, __VA_ARGS__);
#define V(X) V_A(X,NULL)


int
main(int argc, char *argv[])
{
  /* global options */
  int port = 10102;
  int batch_size = 100000;
  int num_batches = -1;
  int nv = 1024;
  int is_int = 0;
  struct hostent * server = NULL;

  int opt = 0;
  while(-1 != (opt = getopt(argc, argv, "p:b:a:x:y:n:i"))) {
    switch(opt) {
      case 'p': {
	port = atoi(optarg);
      } break;

      case 'x': {
	batch_size = atol(optarg);
      } break;

      case 'y': {
	num_batches = atol(optarg);
      } break;

      case 'i': {
	is_int = 1;
      } break;

      case 'a': {
	server = gethostbyname(optarg);
	if(NULL == server) {
	  E_A("ERROR: server %s could not be resolved.", optarg);
	  exit(-1);
	}
      } break;

      case 'n': {
	nv = atol(optarg);
      } break;

      case '?':
      case 'h': {
	printf("Usage:    %s [-p port] [-a server_addr] [-n num_vertices] [-x batch_size] [-y num_batches] [-i]\n", argv[0]);
	printf("Defaults:\n\tport: %d\n\tserver: localhost\n\tnum_vertices: %d\n -i forces the use of integers in place of strings\n", port, nv);
	exit(0);
      } break;
    }
  }

  V_A("Running with: port: %d\n", port);

  /* connect to localhost if server is unspecified */
  if(NULL == server) {
    server = gethostbyname("localhost");
    if(NULL == server) {
      E_A("ERROR: server %s could not be resolved.", "localhost");
      exit(-1);
    }
  }

  /* start the connection */
  int sock_handle = connect_to_batch_server (server, port);

  /* actually generate and send the batches */
  char * buf = NULL, ** fields = NULL;
  uint64_t bufSize = 0, * lengths = NULL, fieldsSize = 0, count = 0;
  int64_t line = 0;
  int batch_num = 0;

  srand (time(NULL));

  while(1) {
    StingerBatch batch;
    batch.set_make_undirected(true);
    batch.set_type(is_int ? NUMBERS_ONLY : STRINGS_ONLY);
    batch.set_keep_alive(true);

    std::string src, dest;

    // IF YOU MAKE THIS LOOP PARALLEL, 
    // MOVE THE SRC/DEST STRINGS ABOVE
    for(int e = 0; e < batch_size; e++) {
      line++;

      int64_t u = rand() % nv;
      int64_t v = rand() % nv;

      if(u == v) {
	e--;
	line--;
	continue;
      }

      /* is insert? */
      EdgeInsertion * insertion = batch.add_insertions();
      if(is_int) {
	insertion->set_source(u);
	insertion->set_destination(v);
      } else {
	insertion->set_source_str(build_name(src, u));
	insertion->set_destination_str(build_name(dest, v));
      }
      insertion->set_weight(1);
      insertion->set_time(line);
    }

    V("Sending messages.");

    send_message(sock_handle, batch);
    sleep(2);

    batch_num++;
    if((batch_num >= num_batches) && (num_batches != -1)) {
      break;
    }
  }

  StingerBatch batch;
  batch.set_make_undirected(true);
  batch.set_type(is_int ? NUMBERS_ONLY : STRINGS_ONLY);
  batch.set_keep_alive(false);
  send_message(sock_handle, batch);

  return 0;
}

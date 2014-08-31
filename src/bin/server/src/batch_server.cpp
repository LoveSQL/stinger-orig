#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <errno.h>

#include "server.h"
#include "stinger_net/stinger_server_state.h"

typedef struct {
  struct stinger * S;
  int sock;
} handle_stream_args;

void *
handle_stream(void * args)
{
  StingerServerState & server_state = StingerServerState::get_server_state();

  struct stinger *S = ((handle_stream_args *) args)->S;
  int sock = ((handle_stream_args *) args)->sock;
  free(args);

  int nfail = 0;

  V("Ready to accept messages.");
  while(1)
  {
    StingerBatch * batch = new StingerBatch();
    if (recv_message(sock, *batch)) {
      nfail = 0;

      V_A("Received message of size %ld", (long)batch->ByteSize());

      if (0 == batch->insertions_size () && 0 == batch->deletions_size ()) {
	V("Empty batch.");
	if (!batch->keep_alive ()) {
	  delete batch;
	  break;
	} else {
	  delete batch;
	  continue;
	}
      }

      server_state.enqueue_batch(batch);

      if(!batch->keep_alive()) {
	delete batch;
	break;
      }

    } else {
      ++nfail;
      V("ERROR Parsing failed.\n");
      if (nfail > 2) break;
    }
  }
}

void
start_tcp_batch_server (struct stinger * S, char * stinger_loc, int port_streams, int port_algs)
{
  StingerServerState & server_state = StingerServerState::get_server_state();
  server_state.set_stinger(S);
  server_state.set_stinger_loc(stinger_loc);
  server_state.set_port(port_algs);
  server_state.set_mon_stinger(stinger_loc, sizeof(stinger_t) + S->length);

  int sock_handle, newsockfd;
  pthread_t garbage_thread_handle;
  socklen_t clilen;
  struct sockaddr_in serv_addr, cli_addr;
  pid_t pid;

  if (-1 == (sock_handle = socket(AF_INET, SOCK_STREAM, 0))) {
    perror("Socket create failed.\n");
    exit(-1);
  }

  bzero ((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(port_streams);

  if (-1 == bind(sock_handle, (struct sockaddr *) &serv_addr, sizeof(serv_addr))) {
    perror("Socket bind failed.\n");
    exit(-1);
  }

  if (-1 == listen(sock_handle, 5)) {
    perror("Socket listen failed.\n");
    exit(-1);
  }

  clilen = sizeof(cli_addr);

  V_A("STINGER server listening for input on port %d",
      (int)port_streams);

  pthread_t alg_handling;
  pthread_create(&alg_handling, NULL, start_alg_handling, NULL);

  while (1) {
    newsockfd = accept (sock_handle, (struct sockaddr *) &cli_addr, &clilen);

    if (newsockfd < 0) {
      perror("Accept new connection failed.\n");
      exit(-1);
    }

    handle_stream_args * args = (handle_stream_args *)calloc(1, sizeof(handle_stream_args));
    args->S = S;
    args->sock = newsockfd;

    pthread_create(&garbage_thread_handle, NULL, handle_stream, (void *)args);
  }

  close(sock_handle);
  return;
}

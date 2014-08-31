#ifndef  STINGER_SHARED_H
#define  STINGER_SHARED_H

#include "stinger.h"

struct stinger_shared;

void *
shmmap (char * name, int oflags, mode_t mode, int prot, size_t size);

int
shmunmap (char * name, void * ptr, size_t size);

struct stinger *
stinger_shared_new (struct stinger_shared ** shared, char ** name);

struct stinger *
stinger_shared_map (struct  stinger_shared ** shared, char * name);

struct stinger *
stinger_shared_free (struct stinger *S, struct stinger_shared * shared, char * name);

struct stinger *
stinger_shared_unmap (struct stinger *S, struct stinger_shared * shared, char * name);

#endif  /*STINGER_SHARED_H*/

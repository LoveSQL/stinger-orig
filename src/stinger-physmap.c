#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "stinger-physmap.h"
#include "stinger-atomics.h"

#define MARKERINT INT64_MAX
uint64_t MARKER_;
void * MARKER = (void *)&MARKER_;

#define CHILDREN_COUNT 256

/* struct defs */

typedef struct tree_node {
  struct tree_node * children[CHILDREN_COUNT];
  struct tree_node * parent;
  uint8_t isEndpoint;
  uint64_t vertexID;
  uint64_t depth;
  char value;
} tree_node_t;

struct stinger_physmap {
  tree_node_t keyTree[MAX_NODES];
  uint64_t keyTreeTop;
  tree_node_t * vtxStack[MAX_VTXID];
  uint64_t vtxStackTop;
};

/* internal functions protos */

tree_node_t *
allocateTreeNode (stinger_physmap_t * map, tree_node_t * parent, uint64_t depth, char value);

int
insertIntoTree(stinger_physmap_t * map, tree_node_t ** node, char * string, uint64_t length);

/* function defs */

/** @brief Allocate and initialize a new physical mapper.
 *
 *  The user is responsible for freeing via stinger_physmap_delete().
 *
 *  @return A new physical mapper.
 */
stinger_physmap_t * 
stinger_physmap_create() {
  stinger_physmap_t * out = malloc(sizeof(stinger_physmap_t));
  if(out) {
    out->vtxStackTop = 0;
    out->keyTreeTop = 0;
    allocateTreeNode(out, NULL, 0, '\0');
  }
  return out;
}

/** @brief Free a physical mapper.
 *
 *  @param map The physical mapper to be freed.
 */
void
stinger_physmap_delete(stinger_physmap_t * map) {
  free(map);
}

tree_node_t *
allocateTreeNode (stinger_physmap_t * map, tree_node_t * parent, uint64_t depth, char value) {
  uint64_t myNode = stinger_int64_fetch_add(&(map->keyTreeTop),1);
  if(map->keyTreeTop >= MAX_NODES) {
    fprintf(stderr, "PHYSMAP: ERROR Out of treenodes\n");
    return NULL;
  }
  memset(map->keyTree + myNode, 0, sizeof(tree_node_t));
  map->keyTree[myNode].parent = parent;
  map->keyTree[myNode].depth = depth;
  map->keyTree[myNode].value = value;
  return map->keyTree + myNode;
}

int
insertIntoTree(stinger_physmap_t * map, tree_node_t ** node, char * string, uint64_t length) {
  if(length == 0) {
    if((*node)->isEndpoint)
      return 1;
    (*node)->isEndpoint = 1;
    if(!stinger_int64_cas(&((*node)->vertexID), 0, MARKERINT)) {
      return 0;
    } else {
      return 2;
    }
  } else {
    if(!(*node)->children[string[0]]) {
      if(!stinger_ptr_cas((void **)&((*node)->children[string[0]]), NULL, MARKER)) {
	(*node)->children[string[0]] = allocateTreeNode(map, *node, (*node)->depth+1, string[0]);
      }
    }
    if(!(*node)->children[string[0]]) {
      return -1;
    }
    while((*node)->children[string[0]] == MARKER) ;
    (*node) = (*node)->children[string[0]];
    return insertIntoTree(map, node, ++string, --length);
  }
}

/** @brief Create a new mapping from a binary data string to a vertex ID.
 *
 *  This function will uniquely map an arbitrary binary string or character
 *  string to a vertex ID in the space of 0 to NV where NV is the number
 *  of unique strings that have been mapped thus far (in other words the vertex
 *  ID space is compact).  It will return -1 on error or if the mapping already exists.
 *  It is safe to call this function in parallel with any other physical mapper function.
 *  To determine if a -1 result is from an error, call stinger_physmap_get_mapping() 
 *  on the same string.  If it also returns -1, then an error has occurred.
 *
 *  @param map The physical mapper.
 *  @param string The binary or character data string.
 *  @param length The length of the string.
 *  @return A unique vertex ID or -1 if the mapping exists or an error occurs.
 */
uint64_t
stinger_physmap_create_mapping (stinger_physmap_t * map, char * string, uint64_t length) {
  if(map->vtxStackTop == MAX_VTXID) {
    fprintf(stderr, "PHYSMAP: ERROR Out of vertices\n");
    return -1;
  }
  uint64_t vertexID;
  tree_node_t * node = map->keyTree;
  int result = insertIntoTree(map, &node, string, length);
  switch(result) {
    case 0:
      vertexID = stinger_int64_fetch_add(&(map->vtxStackTop), 1);
      node->vertexID = vertexID;
      map->vtxStack[vertexID] = node;
      break;
    case 2:
      while(node->vertexID == MARKERINT) ;
      vertexID = node->vertexID;
    default:
      return -1;
  }
  return vertexID;
}

/** @brief Lookup a mapping from a binary data string to a vertex ID.
 *
 *  This function will lookup and return a previously created mapping.  It will return -1 
 *  if no mapping exists.
 *  It is safe to call this function in parallel with any other physical mapper function.
 *
 *  @param map The physical mapper.
 *  @param string The binary or character data string.
 *  @param length The length of the string.
 *  @return A unique vertex ID or -1 if the mapping does not exist.
 */
uint64_t
stinger_physmap_get_mapping (stinger_physmap_t * map, char * string, uint64_t length) {
  tree_node_t * cur = map->keyTree;
  while(length > 0 && cur) {
    cur = cur->children[string[0]];
    string++;
    length--;
  }
  if(cur && cur->isEndpoint) {
    return cur->vertexID;
  } else {
    return -1;
  }
}

/** @brief Lookup the string mapped to a particular vertex ID.
 *
 *  This function will lookup and return a previously created mapping.  It will return -1 
 *  no mapping exists or a reallocation of the output buffer fails.  If the output buffer
 *  is not long enough, this function will reallocate the buffer and update the output buffer
 *  length.
 *  It is safe to call this function in parallel with any other physical mapper function.
 *
 *  @param map The physical mapper.
 *  @param outbuffer A buffer to store the output string.
 *  @param outbufferlength The length of the buffer.
 *  @param vertexID The vertex ID to reverse lookup.
 *  @return 0 on success, -1 on failure.
 */
int
stinger_physmap_get_key (stinger_physmap_t * map, char ** outbuffer, uint64_t * outbufferlength, uint64_t vertexID) {
  if(vertexID >= map->vtxStackTop || map->vtxStack[vertexID] == NULL || (!map->vtxStack[vertexID]->isEndpoint)) 
    return -1;
  tree_node_t * node = map->vtxStack[vertexID];
  if(node->depth+1 > (*outbufferlength)) {
    char * tmpbuffer = realloc(*outbuffer, sizeof(char) * (node->depth+1));
    if(tmpbuffer) {
      (*outbuffer) = tmpbuffer;
      (*outbufferlength) = node->depth+1;
    } else {
      return -1;
    }
  }
  (*outbuffer)[node->depth] = '\0';
  while(node->parent) {
    (*outbuffer)[node->depth-1] = node->value;
    node = node->parent;
  }
  return 0;
}

uint64_t
stinger_physmap_remove_mapping (stinger_physmap_t * map, uint64_t vertexID) {
  printf("***TODO***\n");
}

/* Independent test-driver code */
#if PHYSMAP_TEST

#include <omp.h>
#include <sys/time.h>


struct timeval tv;

double firsttic = 0;
double lasttic = 0;

void tic_reset() {
	gettimeofday(&tv, NULL);
	firsttic = (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
	lasttic = firsttic;
}
double tic_total() {
	gettimeofday(&tv, NULL);
	lasttic = (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
	return lasttic - firsttic;
}
double tic_sincelast() {
	gettimeofday(&tv, NULL);
	double rtnval = ((double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec) - lasttic;
	lasttic = (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
	return rtnval;
}

/*
 * Parallel test driver code
 */

int main(int argc, char *argv[]) {
  stinger_physmap_t * map = stinger_physmap_create();
  if(!map) {
    printf("ALLOC FAILED");
    return 0;
  }

  if(argc < 2) {
    return -1;
  }
  int threads = atoi(argv[1]);
  omp_set_num_threads(threads);
  uint64_t lines_in_file = 0;
  char ** strings;
  uint64_t * lengths;
  float insertion, lookup, reverselookup;
#pragma omp parallel
  {
#pragma omp master
    {
      printf("%d,", omp_get_num_threads());
      FILE * fp = fopen(argv[2], "r");
      char * string = malloc(100*sizeof(char));;
      uint64_t read = 0;
      int bytes = 100;
      while((read = getline(&string, &bytes, fp)) != EOF ) {
	lines_in_file++;
      }
      free(string);
      fclose(fp);
      fp = fopen(argv[2], "r");
      strings = malloc(lines_in_file * sizeof(char *));
      lengths = malloc(lines_in_file * sizeof(uint64_t));
      for(uint64_t i = 0; i < lines_in_file; ++i) {
	string = malloc(100*sizeof(char));;
	read = getline(&string, &bytes, fp);
	strings[i] = string;
	lengths[i] = read - 2;
      }
      printf("%d,",lines_in_file);
    }
  }

  tic_reset();
#pragma omp parallel for
  for(uint64_t i = 0; i < lines_in_file; ++i) {
    uint64_t mapping = stinger_physmap_create_mapping(map, strings[i ], lengths[i ]);
  }
  insertion = tic_sincelast();

#pragma omp parallel for
  for(uint64_t i = 0; i < lines_in_file; ++i) {
    uint64_t mapping = stinger_physmap_get_mapping(map, strings[i ], lengths[i ]);
    if(mapping == -1)
      printf("lu %s %lu %lu\n", strings[i ], lengths[i ], mapping);
  }
  lookup = tic_sincelast();

#pragma omp parallel
{
  char * string2 = malloc(sizeof(char) * 100);
  tic_reset();
#pragma omp for
  for(uint64_t i = 0; i < map->vtxStackTop; ++i) {
    uint64_t slen2 = 100;
    if(stinger_physmap_get_key(map, &string2, &slen2, i ))
      printf("rlu %s %lu\n", string2, slen2);
  }
}
  reverselookup = tic_sincelast();

  printf("%f,%f,%f,", insertion, lookup, reverselookup);
  printf("%f,%f,%f,\n", lines_in_file/insertion, lines_in_file/lookup, lines_in_file/reverselookup);
  stinger_physmap_delete(map);
}

#endif

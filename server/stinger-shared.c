#include "stinger-shared.h"
#include "x86-full-empty.h"
#include "xmalloc.h"

#include <unistd.h>
#include <sys/mman.h>

#define MAX_NAME_LEN 256

struct stinger_shared {
  char ebpool[MAX_NAME_LEN];
  char LVA[MAX_NAME_LEN];
  char ETA[MAX_NAME_LEN];
};

static char * ebpool_name;

/** @brief Wrapper function to open and map shared memory.  
 * 
 * Internally calls shm_open, ftruncate, and mmap.  Can be used to open/map
 * existing shared memory or new shared memory. mmap() flags are always just
 * MAP_SHARED. Ftruncate is silently ignored (for secondary mappings). Names
 * must be less than 256 characters.
 *
 * @param name The string (beginning with /) name of the shared memory object.
 * @param oflags The open flags for shm_open.
 * @param mode The mode flags for shm_open.
 * @param prot The protection/permission flags for mmap. 
 * @param size The size of the memory object to be mapped.
 * @return A pointer to the shared memory or NULL on failure.
 */
void *
shmmap (char * name, int oflags, mode_t mode, int prot, size_t size) 
{
#if !defined(__MTA__)
  int fd = shm_open(name, oflags, mode);
#else
  int fd = open(name, oflags);
#endif

  if(fd == -1) {
    fprintf(stderr, "\nSHMMAP shm_open ERROR %s\n", strerror(errno)); fflush(stdout);
    return NULL;
  } 

#if !defined(__MTA__)
  int dontcare = ftruncate(fd, size);
  /* silently ignore ftruncate errors */
  
  void * rtn = mmap(NULL, size, prot, MAP_SHARED, fd, 0);
#else
  void * rtn = mmap(NULL, size, prot, MAP_SHARED|MAP_ANON, fd, 0);
#endif

  if(rtn == MAP_FAILED) {
    fprintf(stderr, "\nSHMMAP mmap ERROR %s\n", strerror(errno)); fflush(stdout);
    return NULL;
  }

  return rtn;
} 

/** @brief Wrapper function to unmap / unlink shared memory.
 * 
 * @param name The string (beginning with /) name of the shared memory object.
 * @param ptr The pointer to the shared memory.
 * @param size The size of the memory object to be removed.
 * @return 0 on success, -1 on failure.
 */
int
shmunmap (char * name, void * ptr, size_t size) 
{
  if(munmap(ptr, size))
    return -1;

#if !defined(__MTA__)
  if(shm_unlink(name))
    return -1;
#else
  if(unlink(name))
    return -1;
#endif

  return 0;
}

/** @brief Initialize the ebpool as a shared memory object.
 */
static void 
init_shared_ebpool (void) 
{
  struct stinger_eb *new_ebpool = (struct stinger_eb *)readfe ((uint64_t *)&ebpool);
  if (new_ebpool) {
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
    return;
  } else {
    ebpool_name = xmalloc(MAX_NAME_LEN * sizeof(char));
#if !defined(__MTA__)
    sprintf(ebpool_name, "/%lx", (uint64_t)rand());
#else
    char *pwd = xmalloc (sizeof(char) * (MAX_NAME_LEN-16));
    getcwd(pwd,  MAX_NAME_LEN-16);
    sprintf(ebpool_name, "%s/%lx", pwd, (uint64_t)rand());
    free(pwd);
#endif
    new_ebpool = shmmap (ebpool_name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
      PROT_READ | PROT_WRITE, EBPOOL_SIZE * sizeof (struct stinger_eb));
    ebpool_tail = 1;
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
  }
}

/** @brief Map in a shared ebpool.
 *
 * @param shared The stinger_shared struct containing the mapping name.
 */
static void 
map_shared_ebpool (struct stinger_shared * shared) 
{
  struct stinger_eb *new_ebpool = (struct stinger_eb *)readfe ((uint64_t *)&ebpool);
  if (new_ebpool) {
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
    return;
  } else {
    ebpool_name = xmalloc(MAX_NAME_LEN * sizeof(char));
    memcpy(ebpool_name, shared->ebpool, MAX_NAME_LEN);
    new_ebpool = shmmap (ebpool_name, O_RDONLY, S_IRUSR | S_IWUSR, 
      PROT_READ, EBPOOL_SIZE * sizeof (struct stinger_eb));
    ebpool_tail = 1;
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
  }
}

/** @brief Unmap and unlink the shared ebpool.
 */
static void
free_shared_ebpool (void)
{
  shmunmap (ebpool_name, ebpool, EBPOOL_SIZE);
  ebpool = NULL;
  ebpool_tail = EBPOOL_SIZE; /* to prevent getting eb's from empty pool */
}

/** @brief Create a new empty STINGER in shared memory.
 *
 * The stinger_shared is a struct containing all of the strings needed to map
 * in the resulting STINGER from another program.  The shared STINGER is stored
 * in shared memory mapped using the string stored in out.  The out string is 
 * all that is needed to map in the same stinger in another program.
 *
 * @param S The return for the stinger_shared structure.
 * @param out The return for the random stinger_shared name.
 * @return A pointer to the new stinger.
 */
struct stinger *
stinger_shared_new (struct stinger_shared ** S, char ** out)
{
  *out = xmalloc(sizeof(char) * MAX_NAME_LEN);
#if !defined(__MTA__)
  sprintf(*out, "/%lx", (uint64_t)rand());
#else
  char *pwd = xmalloc (sizeof(char) * (MAX_NAME_LEN-16));
  getcwd(pwd,  MAX_NAME_LEN-16);
  sprintf(*out, "%s/%lx", pwd, (uint64_t)rand());
#endif
  struct stinger_shared * shared = shmmap(*out, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR,
    PROT_READ | PROT_WRITE, 1 * sizeof(*shared));
  struct stinger *G = xcalloc (1, sizeof (*G));
  size_t i;

  if (!readff((uint64_t *)&ebpool))
    init_shared_ebpool();

  memcpy(shared->ebpool, ebpool_name, MAX_NAME_LEN);
  printf("\n\tshared->ebpool %s", shared->ebpool);

#if !defined(__MTA__)
  sprintf(shared->LVA, "/%lx", (uint64_t)rand());
#else
  sprintf(shared->LVA, "%s/%lx", pwd, (uint64_t)rand());
#endif
  G->LVA = shmmap(shared->LVA, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
    PROT_READ | PROT_WRITE, STINGER_MAX_LVERTICES * sizeof (G->LVA[0]));
  xzero(G->LVA, STINGER_MAX_LVERTICES * sizeof(G->LVA[0]));
  G->LVASize = STINGER_MAX_LVERTICES;

#if !defined(__MTA__)
  sprintf(shared->ETA, "/%lx", (uint64_t)rand());
#else
  sprintf(shared->ETA, "%s/%lx", pwd, (uint64_t)rand());
  free(pwd);
#endif
  G->ETA = shmmap(shared->ETA, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
    PROT_READ | PROT_WRITE, STINGER_NUMETYPES * sizeof(struct stinger_etype_array));

  //OMP ("omp parallel for")
  MTA ("mta assert parallel")
  MTASTREAMS ()
  for (i = 0; i < STINGER_NUMETYPES; ++i) {
    G->ETA[i].high = 0;
    G->ETA[i].length = EBPOOL_SIZE;
  }

  *S = shared;
  return G;
}

/** @brief Map in an existinger STINGER in shared memory.
 *
 * Input the name obtained from sting_shared_new above in a different program
 * to map in the same STINGER in shared memory.
 *
 * @param S The return for the stinger_shared structure.
 * @param out The input for the stinger_shared name.
 * @return A pointer to the stinger.
 */
struct stinger *
stinger_shared_map (struct  stinger_shared ** shared, char * name) {
  *shared = shmmap(name, O_RDONLY, S_IRUSR,
    PROT_READ, 1 * sizeof(struct stinger_shared));
  struct stinger *G = xcalloc (1, sizeof (*G));

  if (!readff((uint64_t *)&ebpool))
    map_shared_ebpool(*shared);

  G->LVA = shmmap((*shared)->LVA, O_RDONLY, S_IRUSR, 
    PROT_READ, STINGER_MAX_LVERTICES * sizeof (G->LVA[0]));
  G->LVASize = STINGER_MAX_LVERTICES;

  G->ETA = shmmap((*shared)->ETA, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR, 
    PROT_READ | PROT_WRITE, STINGER_NUMETYPES * sizeof(struct stinger_etype_array));

  return G;
}

/** @brief Unmap and unlink a shared STINGER from stinger_shared_new.
 * 
 * @param S The STINGER pointer.
 * @param shared The stinger_shared pointer.
 * @param name The name returned from stinger_shared_new.
 * @return A NULL pointer. Name and shared will also be freed.
 */
struct stinger *
stinger_shared_free (struct stinger *S, struct stinger_shared * shared, char * name)
{
  if (!S)
    return S;

  shmunmap (shared->ETA, S->ETA, STINGER_NUMETYPES * sizeof(struct stinger_etype_array));
  shmunmap (shared->LVA, S->LVA, STINGER_MAX_LVERTICES * sizeof(S->LVA[0]));
  free (S);
  shmunmap(name, shared, sizeof(*shared));
  free_shared_ebpool ();
  free (name);
  return NULL;
}

/** @brief Unmap a shared STINGER from another program.
 * 
 * @param S The STINGER pointer.
 * @param shared The stinger_shared pointer.
 * @param name The name used to map the shared STINGER.
 * @return A NULL pointer.
 */
struct stinger *
stinger_shared_unmap (struct stinger *S, struct stinger_shared * shared, char * name)
{
  if(!S)
    return S;

  /* 
   * Letting the program closing clean these up since it seems to cause 
   * problems otherwise 
   */

  free(S);
  return NULL;
}

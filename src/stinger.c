/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#include <dirent.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "stinger.h"
#include "stinger-atomics.h"
#include "stinger-utils.h"
#include "xmalloc.h"
#include "x86-full-empty.h"

/* {{{ Edge block pool */
/* TODO XXX Rework / possibly move EB POOL functions */

uint64_t ebpool_tail = 0;
struct stinger_eb *ebpool = NULL;

MTA ("mta inline")
static void init_ebpool (void)
{
  struct stinger_eb *new_ebpool = (struct stinger_eb *)readfe ((uint64_t *)&ebpool);
  if (new_ebpool) {
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
    return;
  } else {
    new_ebpool = xmalloc (EBPOOL_SIZE * sizeof (struct stinger_eb));
    ebpool_tail = 1;
    writeef ((uint64_t *)&ebpool, (uint64_t)new_ebpool);
  }
}

static void
free_ebpool (void)
{
  free (ebpool);
  ebpool = NULL;
  ebpool_tail = EBPOOL_SIZE; /* to prevent getting eb's from empty pool */
}

MTA ("mta expect parallel context") MTA ("mta inline")
static void
get_from_ebpool (eb_index_t *out, size_t k)
{
  eb_index_t ebt0;
  {
    ebt0 = stinger_size_fetch_add (&ebpool_tail, k);
    if (ebt0 + k >= EBPOOL_SIZE) {
      fprintf (stderr, "XXX: eb pool exhausted\n");
      abort ();
    }
    OMP("omp parallel for")
      MTA ("mta assert nodep")
      MTASTREAMS ()MTA ("mta block schedule")
      for (size_t ki = 0; ki < k; ++ki)
        out[ki] = ebt0 + ki;
  }
}

/* }}} */

/* {{{ Internal utilities */

/** @brief Returns the out-degree of a vertex
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @return Out-degree of vertex i_
 */
int64_t
stinger_outdegree (const struct stinger *S_, int64_t i_)
{
  return S_->LVA[i_].outDegree;
}

/** @brief Returns the out-degree of a vertex for a given edge type
 *
 *  @param S The STINGER data structure
 *  @param i Logical vertex ID
 *  @param type Edge type
 *  @return Out-degree of vertex i with type
 */
int64_t
stinger_typed_outdegree (const struct stinger * S, int64_t i, int64_t type) {
  int64_t out = 0;
  struct stinger_eb *current_eb = ebpool + (S)->LVA[(i)].edges;
  while(current_eb != ebpool) {
    if(current_eb->etype == type) {
      out += current_eb->numEdges;
    }
    current_eb = ebpool + current_eb->next;
  }
  return out;
}

/** @brief Calculate the largest active vertex ID
 *
 *  Finds the largest vertex ID whose in-degree and/or out-degree
 *  is greater than zero.
 *
 *  <em>NOTE:</em> If you are using this to obtain a
 *  value of "nv" for additional STINGER calls, you must add one to the 
 *  result.
 *
 *  @param S The STINGER data structure
 *  @return Largest active vertex ID
 */
uint64_t
stinger_max_active_vertex(const struct stinger * S) {
  uint64_t out = 0;
  OMP("omp parallel") {
    uint64_t local_max = 0;
    OMP("omp for")
    for(uint64_t i = 0; i < STINGER_MAX_LVERTICES; i++) {
      if((stinger_indegree(S, i) > 0 || stinger_outdegree(S, i) > 0) && 
	i > local_max) {
	local_max = i;
      }
    }
    OMP("omp critical") {
      if(local_max > out)
	out = local_max;
    }
  }
  return out;
}

/** @brief Calculate the number of active vertices
 *
 *  Counts the number of vertices whose in-degree is greater than zero and
 *  out-degree is greater than zero.
 *
 *  @param S The STINGER data structure
 *  @return Number of active vertices
 */
uint64_t
stinger_num_active_vertices(const struct stinger * S) {
  uint64_t out = 0;
  OMP("omp parallel for reduction(+:out)")
  for(uint64_t i = 0; i < STINGER_MAX_LVERTICES; i++) {
    if(stinger_indegree(S, i) > 0 || stinger_outdegree(S, i) > 0) {
      out++;
    }
  }
  return out;
}

/** @brief Returns the in-degree of a vertex
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @return In-degree of vertex i_
 */
int64_t
stinger_indegree (const struct stinger * S_, int64_t i_)
{
  return S_->LVA[i_].inDegree;
}

/** @brief Returns the vertex weight
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @return Weight of vertex i_
 */
int64_t
stinger_vweight (const struct stinger * S_, int64_t i_)
{
  return S_->LVA[i_].weight;
}

/** @brief Sets the vertex weight
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @param weight_ Vertex weight
 *  @return 1 on success
 */
int64_t
stinger_set_vweight (struct stinger * S_, int64_t i_, int64_t weight_)
{
  return (S_->LVA[i_].weight = weight_);
}

/** @brief Returns the vertex type
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @return Type of vertex i_
 */
int64_t
stinger_vtype (const struct stinger * S_, int64_t i_)
{
  return S_->LVA[i_].vtype;
}

/** @brief Sets the vertex type
 *
 *  @param S_ The STINGER data structure
 *  @param i_ Logical vertex ID
 *  @param type_ Vertex type
 *  @return 1 on success
 */
int64_t
stinger_set_vtype (const struct stinger * S_, int64_t i_, int64_t type_)
{
  return (S_->LVA[i_].vtype = type_);
}

const struct stinger_eb *
stinger_edgeblocks (const struct stinger *S_, int64_t i_)
{
  return ebpool + S_->LVA[i_].edges;
}

const struct stinger_eb *
stinger_next_eb (const struct stinger *G /*UNUSED*/,
                 const struct stinger_eb *eb_)
{
  return ebpool + readff((uint64_t *)&eb_->next);
}

int64_t
stinger_eb_type (const struct stinger_eb * eb_)
{
  return eb_->etype;
}

int
stinger_eb_high (const struct stinger_eb *eb_)
{
  return eb_->high;
}

int
stinger_eb_is_blank (const struct stinger_eb *eb_, int k_)
{
  return eb_->edges[k_].neighbor < 0;
}

int64_t
stinger_eb_adjvtx (const struct stinger_eb * eb_, int k_)
{
  return eb_->edges[k_].neighbor;
}

int64_t
stinger_eb_weight (const struct stinger_eb * eb_, int k_)
{
  return eb_->edges[k_].weight;
}

int64_t
stinger_eb_ts (const struct stinger_eb * eb_, int k_)
{
  return eb_->edges[k_].timeRecent;
}

int64_t
stinger_eb_first_ts (const struct stinger_eb * eb_, int k_)
{
  return eb_->edges[k_].timeFirst;
}

/**
* @brief Count the total number of edges in STINGER.
*
* @param S The STINGER data structure
*
* @return The number of edges in STINGER
*/
int64_t
stinger_total_edges (const struct stinger * S)
{
  uint64_t rtn = 0;
  for (uint64_t i = 0; i < STINGER_MAX_LVERTICES; i++) {
    rtn += stinger_outdegree (S, i);
  }
  return rtn;
}

/**
* @brief Calculate the total size of the active STINGER graph in memory.
*
* @param S The STINGER data structure
*
* @return The number of bytes currently in use by STINGER
*/
size_t
stinger_graph_size (const struct stinger *S)
{
  int64_t num_edgeblocks = S->ETA->high;
  int64_t size_edgeblock = sizeof(struct stinger_eb);

  int64_t size_vb = sizeof(struct stinger_vb);
  int64_t num_vertexblocks = S->LVASize;

  int64_t result = (num_edgeblocks * size_edgeblock) + (num_vertexblocks * size_vb);

  return result;
}

void
stinger_print_eb(struct stinger_eb * eb) {
  printf(
    "EB VTX:  %ld\n"
    "  NEXT:    %ld\n"
    "  TYPE:    %ld\n"
    "  NUMEDGS: %ld\n"
    "  HIGH:    %ld\n"
    "  SMTS:    %ld\n"
    "  LGTS:    %ld\n"
    "  EDGES:\n",
    eb->vertexID, eb->next, eb->etype, eb->numEdges, eb->high, eb->smallStamp, eb->largeStamp);
  uint64_t j = 0;
  for (; j < eb->high && j < STINGER_EDGEBLOCKSIZE; j++) {
    printf("    TO: %s%ld WGT: %ld TSF: %ld TSR: %ld\n", 
      eb->edges[j].neighbor < 0 ? "x " : "  ", eb->edges[j].neighbor < 0 ? ~(eb->edges[j].neighbor) : eb->edges[j].neighbor, 
      eb->edges[j].weight, eb->edges[j].timeFirst, eb->edges[j].timeRecent);
  }
  if(j < STINGER_EDGEBLOCKSIZE) {
    printf("  ABOVE HIGH:\n");
    for (; j < STINGER_EDGEBLOCKSIZE; j++) {
      printf("    TO: %ld WGT: %ld TSF: %ld TSR: %ld\n", 
	eb->edges[j].neighbor, eb->edges[j].weight, eb->edges[j].timeFirst, eb->edges[j].timeRecent);
    }
  }
}

/**
* @brief Checks the STINGER metadata for inconsistencies.
*
* @param S The STINGER data structure
* @param NV The total number of vertices
*
* @return 0 on success; failure otherwise
*/
uint32_t
stinger_consistency_check (struct stinger *S, uint64_t NV)
{
  uint32_t returnCode = 0;
  uint64_t *inDegree = calloc (NV, sizeof (uint64_t));
  if (inDegree == NULL) {
    returnCode |= 0x00000001;
    return returnCode;
  }

  // check blocks
  OMP("omp parallel for reduction(|:returnCode)")
  MTA("mta assert nodep")
  for (uint64_t i = 0; i < NV; i++) {
    uint64_t curOutDegree = 0;
    const struct stinger_eb *curBlock = stinger_edgeblocks(S, i);
    while (curBlock != ebpool) {
      if (curBlock->vertexID != i)
        returnCode |= 0x00000002;
      if (curBlock->high > STINGER_EDGEBLOCKSIZE)
        returnCode |= 0x00000004;

      int64_t numEdges = 0;
      int64_t smallStamp = INT64_MAX;
      int64_t largeStamp = INT64_MIN;

      uint64_t j = 0;
      for (; j < curBlock->high && j < STINGER_EDGEBLOCKSIZE; j++) {
        if (!stinger_eb_is_blank (curBlock, j)) {
          stinger_int64_fetch_add (&inDegree[stinger_eb_adjvtx (curBlock, j)], 1);
          curOutDegree++;
          numEdges++;
          if (stinger_eb_ts (curBlock, j) < smallStamp)
            smallStamp = stinger_eb_ts (curBlock, j);
          if (stinger_eb_first_ts (curBlock, j) < smallStamp)
            smallStamp = stinger_eb_first_ts (curBlock, j);
          if (stinger_eb_ts (curBlock, j) > largeStamp)
            largeStamp = stinger_eb_ts (curBlock, j);
          if (stinger_eb_first_ts (curBlock, j) > largeStamp)
            largeStamp = stinger_eb_first_ts (curBlock, j);
        }
      }
      if (numEdges && numEdges != curBlock->numEdges)
        returnCode |= 0x00000008;
      if (numEdges && largeStamp > curBlock->largeStamp)
        returnCode |= 0x00000010;
      if (numEdges && smallStamp < curBlock->smallStamp)
        returnCode |= 0x00000020;
      for (; j < STINGER_EDGEBLOCKSIZE; j++) {
        if (!(stinger_eb_is_blank (curBlock, j) ||
              (stinger_eb_adjvtx (curBlock, j) == 0
               && stinger_eb_weight (curBlock, j) == 0
               && stinger_eb_ts (curBlock, j) == 0
               && stinger_eb_first_ts (curBlock, j) == 0)))
          returnCode |= 0x00000040;
      }
      curBlock = ebpool + curBlock->next;
    }

    if (curOutDegree != S->LVA[i].outDegree)
      returnCode |= 0x00000080;
  }

  OMP("omp parallel for reduction(|:returnCode)")
  MTA("mta assert nodep")
  for (uint64_t i = 0; i < NV; i++) {
    if (inDegree[i] != S->LVA[i].inDegree)
      returnCode |= 0x00000100;
  }

  free (inDegree);

#if STINGER_NUMETYPES == 1
  // check for self-edges and duplicate edges
  int64_t count_self = 0;
  int64_t count_duplicate = 0;

  int64_t * off = NULL;
  int64_t * ind = NULL;

  stinger_to_sorted_csr (S, NV, &off, &ind, NULL, NULL, NULL, NULL);

  MTA ("mta assert nodep")
  OMP ("omp parallel for reduction(+:count_self, count_duplicate)")
  for (int64_t k = 0; k < NV; k++)
  {
    int64_t myStart = off[k];
    int64_t myEnd = off[k+1];

    for (int64_t j = myStart; j < myEnd; j++)
    {
      if (ind[j] == k)
	count_self++;
      if (ind[j] == ind[j+1] && j < myEnd-1)
	count_duplicate++;
    }
  }

  free (ind);
  free (off);

  if (count_self != 0)
    returnCode |= 0x00000200;

  if (count_duplicate != 0)
    returnCode |= 0x00000400;
#endif

  return returnCode;
}

/**
* @brief Calculate statistics on edge block fragmentation in the graph
*
* @param S The STINGER data structure
* @param NV The total number of vertices
*
* @return Void
*/
void
stinger_fragmentation (struct stinger *S, uint64_t NV, struct stinger_fragmentation_t * stats)
{
  uint64_t numSpaces = 0;
  uint64_t numBlocks = 0;
  uint64_t numEdges = 0;

  OMP ("omp parallel for reduction(+:numSpaces, numBlocks, numEdges)")
  for (uint64_t i = 0; i < NV; i++) {
    const struct stinger_eb *curBlock = stinger_edgeblocks(S, i);

    while (curBlock != ebpool) {
      uint64_t found = 0;
      for (uint64_t j = 0; j < curBlock->high && j < STINGER_EDGEBLOCKSIZE; j++) {
        if (stinger_eb_is_blank (curBlock, j)) {
          numSpaces++;
	  found = 1;
        }
	else {
	  numEdges++;
	}
      }
      numBlocks += found;
      curBlock = ebpool + curBlock->next;
    }
  }

  int64_t totalEdgeBlocks = S->ETA->high;

  stats->num_empty_edges = numSpaces;
  stats->num_fragmented_blocks = numBlocks;
  stats->num_edges = numEdges;
  stats->edge_blocks_in_use = totalEdgeBlocks;


  double fillPercent = (double) numEdges / (double) (totalEdgeBlocks * STINGER_EDGEBLOCKSIZE);
  stats->fill_percent = fillPercent;
}

/* }}} */



/* {{{ Allocating and tearing down */

MTA ("mta inline")
void stinger_init (void)
{
  if (!ebpool)
    init_ebpool ();
}


/** @brief Create a new STINGER data structure.
 *
 *  Allocates memory for a STINGER data structure.  If this is the first STINGER
 *  to be allocated, it also initializes the edge block pool.  Edge blocks are
 *  allocated and assigned for each value less than STINGER_NUMETYPES.
 *
 *  @return Pointer to struct stinger
 */
MTA ("mta inline")
struct stinger *stinger_new (void)
{
  struct stinger *G = xcalloc (1, sizeof (*G));
  size_t i;

  if (!readff((uint64_t *)&ebpool))
    stinger_init ();

  G->LVA = xcalloc (STINGER_MAX_LVERTICES, sizeof (G->LVA[0]));
  G->LVASize = STINGER_MAX_LVERTICES;
  G->ETA = xmalloc (STINGER_NUMETYPES * sizeof(struct stinger_etype_array));

#if STINGER_NUMETYPES == 1
  G->ETA[0].length = EBPOOL_SIZE;
  G->ETA[0].high = 0;
#else
  OMP ("omp parallel for")
  MTA ("mta assert parallel")
  MTASTREAMS ()
  for (i = 0; i < STINGER_NUMETYPES; ++i) {
    G->ETA[i].length = EBPOOL_SIZE;
    G->ETA[i].high = 0;
  }
#endif

  return G;
}

/** @brief Free memory allocated to a particular STINGER instance.
 *
 *  Frees the ETA pointers for each edge type, the LVA, and the struct stinger
 *  itself.  Does not actually free any edge blocks, as there may be other
 *  active STINGER instances.
 *
 *  @param S The STINGER data structure
 *  @return NULL on success
 */
struct stinger *
stinger_free (struct stinger *S)
{
  size_t i;
  if (!S)
    return S;

  free (S->ETA);
  free (S->LVA);
  free (S);
  free_ebpool ();
  return NULL;
}

/** @brief Free the STINGER data structure and all edge blocks.
 *
 *  Free memory allocated to the specified STINGER data structure.  Also frees
 *  the STINGER edge block pool, effectively ending all STINGER operations globally.
 *  Only call this function if you are done using STINGER entirely.
 *
 *  @param S The STINGER data structure
 *  @return NULL on success
 */
struct stinger *
stinger_free_all (struct stinger *S)
{
  struct stinger *out;
  out = stinger_free (S);
  return out;
}

/* TODO inspect possibly move out with other EB POOL stuff */
MTA ("mta expect parallel context")
static eb_index_t new_eb (int64_t etype, int64_t from)
{
  size_t k;
  eb_index_t out = 0;
  get_from_ebpool (&out, 1);
  struct stinger_eb * block = ebpool + out;
  assert (block != ebpool);
  xzero (block, sizeof (*block));
  block->etype = etype;
  block->vertexID = from;
  block->smallStamp = INT64_MAX;
  block->largeStamp = INT64_MIN;
  return out;
}

MTA ("mta expect parallel context")
void
new_ebs (eb_index_t *out, size_t neb, int64_t etype,
         int64_t from)
{
  if (neb < 1)
    return;
  get_from_ebpool (out, neb);

  OMP ("omp parallel for")
    //MTA("mta assert nodep")
    MTASTREAMS ()MTA ("mta block schedule")
    //MTA("mta parallel single processor")
    for (size_t i = 0; i < neb; ++i) {
      struct stinger_eb * block = ebpool + out[i];
      xzero (block, sizeof (*block));
      block->etype = etype;
      block->vertexID = from;
      block->smallStamp = INT64_MAX;
      block->largeStamp = INT64_MIN;
    }
}

MTA ("mta expect serial context")
static void
new_blk_ebs (eb_index_t *out, const struct stinger *restrict G,
             const int64_t nvtx, const size_t * restrict blkoff,
             const int64_t etype)
{
  size_t neb;
  if (nvtx < 1)
    return;
  neb = blkoff[nvtx];
  get_from_ebpool (out, neb);

  OMP ("omp parallel for")
    MTA ("mta assert nodep")
    MTASTREAMS ()MTA ("mta block schedule")
    for (size_t k = 0; k < neb; ++k) {
      struct stinger_eb * block = ebpool + out[k];
      xzero (block, sizeof (*block));
      block->etype = etype;
      block->smallStamp = INT64_MAX;
      block->largeStamp = INT64_MIN;
    }

  OMP ("omp parallel for")
    MTA ("mta assert nodep")
    MTASTREAMS ()MTA ("mta interleave schedule")
    for (int64_t v = 0; v < nvtx; ++v) {
      const int64_t from = v;
      const size_t blkend = blkoff[v + 1];
      MTA ("mta assert nodep")
        for (size_t k = blkoff[v]; k < blkend; ++k)
          ebpool[out[k]].vertexID = from;
      if (blkend)
        MTA ("mta assert nodep")
          for (size_t k = blkoff[v]; k < blkend - 1; ++k)
            ebpool[out[k]].next = out[k + 1];
    }
}

/* }}} */

/* TODO XXX Clean, remove realloc */
MTA ("mta inline")
void
push_ebs (struct stinger *G, size_t neb,
          eb_index_t * restrict eb)
/* XXX: Reallocing here will nuke current readers *AND* writers.
   A DCAS could save the writers, but we need RCU-like machinery
   for readers. */
{
  int64_t etype, place;
  assert (G);
  assert (eb);

  if (!neb)
    return;

  etype = ebpool[eb[0]].etype;
  assert (etype >= 0);
  assert (etype < STINGER_NUMETYPES);

  place = stinger_int64_fetch_add (&(G->ETA[etype].high), neb);

  eb_index_t *blocks;
  blocks = G->ETA[etype].blocks;

  MTA ("mta assert nodep")
  for (int64_t k = 0; k < neb; ++k)
    blocks[place + k] = eb[k];
}

MTA ("mta inline")
struct curs
etype_begin (struct stinger_vb *v, int etype)
{
  struct curs out;
  assert (v);
  out.eb = readff((uint64_t *)&(v->edges));
  out.loc = &(v->edges);
  while (out.eb && ebpool[out.eb].etype != etype) {
    out.loc = &(ebpool[out.eb].next);
    out.eb = readff((uint64_t *)&(ebpool[out.eb].next));
  }
  return out;
}

MTA ("mta inline")
void
update_edge_data (struct stinger * S, struct stinger_eb *eb,
                  uint64_t index, int64_t neighbor, int64_t weight, int64_t ts)
{
  struct stinger_edge * e = eb->edges + index;

  /* insertion */
  if(neighbor >= 0) {
    e->weight = weight;
    /* is this a new edge */
    if(e->neighbor < 0 || index >= eb->high) {
      e->neighbor = neighbor;

      /* only edge in block? - assuming we have block effectively locked */
      if(stinger_int64_fetch_add(&eb->numEdges, 1) == 0) {
	eb->smallStamp = ts;
	eb->largeStamp = ts;
      }

      /* register new edge */
      stinger_int64_fetch_add (&(S->LVA[eb->vertexID].outDegree), 1);
      stinger_int64_fetch_add (&(S->LVA[neighbor].inDegree), 1);

      if (index >= eb->high)
	eb->high = index + 1;

      writexf(&e->timeFirst, ts);
    }

    /* check metadata and update - lock metadata for safety */
    if (ts < readff(&eb->smallStamp) || ts > eb->largeStamp) {
      int64_t smallStamp = readfe(&eb->smallStamp);
      if (ts < smallStamp)
	smallStamp = ts;
      if (ts > eb->largeStamp)
	eb->largeStamp = ts;
      writeef(&eb->smallStamp, smallStamp);
    }

    e->timeRecent = ts;

  } else if(e->neighbor >= 0) {
    /* are we deleting an edge */
    stinger_int64_fetch_add (&(S->LVA[eb->vertexID].outDegree), -1);
    stinger_int64_fetch_add (&(S->LVA[~neighbor].inDegree), -1);
    stinger_int64_fetch_add (&(eb->numEdges), -1);
    e->neighbor = neighbor;
  } 

  /* we always do this to update weight and  unlock the edge if needed */
}

/** @brief Insert a directed edge.
 *
 *  Inserts a typed, directed edge.  First timestamp is set, if the edge is
 *  new.  Recent timestamp is updated.  Weight is set to specified value regardless.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param weight Edge weight
 *  @param timestamp Edge timestamp
 *  @return 1 if edge is inserted successfully for the first time, 0 if edge is already found and updated, -1 if error.
 */
MTA ("mta inline") MTA("mta serial")
int
stinger_insert_edge (struct stinger *G,
                     int64_t type, int64_t from, int64_t to,
                     int64_t weight, int64_t timestamp)
{
  if(from == to)
    return -1;

  /* Do *NOT* call this concurrently with different edge types. */
  STINGERASSERTS ();

  struct curs curs;
  struct stinger_eb *tmp;
  struct stinger_eb *ebpool_priv = ebpool;

  curs = etype_begin (&G->LVA[from], type);
  /*
  Possibilities:
  1: Edge already exists and only needs updated.
  2: Edge does not exist, fits in an existing block.
  3: Edge does not exist, needs a new block.
  */

  /* 1: Check if the edge already exists. */
  for (tmp = ebpool_priv + curs.eb; tmp != ebpool_priv; tmp = ebpool_priv + readff((uint64_t *)&tmp->next)) {
    if(type == tmp->etype) {
      size_t k, endk;
      endk = tmp->high;

      for (k = 0; k < endk; ++k) {
	if (to == tmp->edges[k].neighbor) {
	  update_edge_data (G, tmp, k, to, weight, timestamp);
	  return 0;
	}
      }
    }
  }

  while (1) {
    eb_index_t * block_ptr = curs.loc;
    curs.eb = readff((uint64_t *)curs.loc);
    /* 2: The edge isn't already there.  Check for an empty slot. */
    for (tmp = ebpool_priv + curs.eb; tmp != ebpool_priv; tmp = ebpool_priv + readff((uint64_t *)&tmp->next)) {
      if(type == tmp->etype) {
	size_t k, endk;
	endk = tmp->high;

	for (k = 0; k < STINGER_EDGEBLOCKSIZE; ++k) {
	  int64_t myNeighbor = tmp->edges[k].neighbor;
	  if (to == myNeighbor && k < endk) {
	    update_edge_data (G, tmp, k, to, weight, timestamp);
	    return 0;
	  }

	  if (myNeighbor < 0 || k >= endk) {
	    int64_t timefirst = readfe ( &(tmp->edges[k].timeFirst) );
	    int64_t thisEdge = tmp->edges[k].neighbor;
	    endk = tmp->high;

	    if (thisEdge < 0 || k >= endk) {
	      update_edge_data (G, tmp, k, to, weight, timestamp);
	      return 1;
	    } else if (to == thisEdge) {
	      update_edge_data (G, tmp, k, to, weight, timestamp);
	      writexf ( &(tmp->edges[k].timeFirst), timefirst);
	      return 0;
	    } else {
	      writexf ( &(tmp->edges[k].timeFirst), timefirst);
	    }
	  }
	}
      }
      block_ptr = &(tmp->next);
    }

    /* 3: Needs a new block to be inserted at end of list. */
    eb_index_t old_eb = readfe ((uint64_t *)block_ptr );
    if (!old_eb) {
      eb_index_t newBlock = new_eb (type, from);
      if (newBlock == 0) {
	writeef ((uint64_t *)block_ptr, (uint64_t)old_eb);
	return -1;
      } else {
	update_edge_data (G, ebpool_priv + newBlock, 0, to, weight, timestamp);
	ebpool_priv[newBlock].next = 0;
	push_ebs (G, 1, &newBlock);
      }
      writeef ((uint64_t *)block_ptr, (uint64_t)newBlock);
      return 1;
    }
    writeef ((uint64_t *)block_ptr, (uint64_t)old_eb);
  }
}

/** @brief Increments a directed edge.
 *
 *  Increments the weight of a typed, directed edge.
 *  Recent timestamp is updated.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param weight Edge weight
 *  @param timestamp Edge timestamp
 *  @return 1
 */
MTA ("mta inline")
int
stinger_incr_edge (struct stinger *G,
                   int64_t type, int64_t from, int64_t to,
                   int64_t weight, int64_t timestamp)
{
  if(from == to)
    return -1;

  /* Do *NOT* call this concurrently with different edge types. */
  STINGERASSERTS ();

  struct curs curs;
  struct stinger_eb *tmp;
  struct stinger_eb *ebpool_priv = ebpool;

  curs = etype_begin (&G->LVA[from], type);
  /*
  Possibilities:
  1: Edge already exists and only needs updated.
  2: Edge does not exist, fits in an existing block.
  3: Edge does not exist, needs a new block.
  */

  /* 1: Check if the edge already exists. */
  for (tmp = ebpool_priv + curs.eb; tmp != ebpool_priv; tmp = ebpool_priv + readff((uint64_t *)&tmp->next)) {
    if(type == tmp->etype) {
      size_t k, endk;
      endk = tmp->high;

      for (k = 0; k < endk; ++k) {
	if (to == tmp->edges[k].neighbor) {
	  update_edge_data (G, tmp, k, to, tmp->edges[k].weight + weight, timestamp);
	  return 0;
	}
      }
    }
  }

  while (1) {
    eb_index_t * block_ptr = curs.loc;
    curs.eb = readff((uint64_t *)curs.loc);
    /* 2: The edge isn't already there.  Check for an empty slot. */
    for (tmp = ebpool_priv + curs.eb; tmp != ebpool_priv; tmp = ebpool_priv + readff((uint64_t *)&tmp->next)) {
      if(type == tmp->etype) {
	size_t k, endk;
	endk = tmp->high;

	for (k = 0; k < STINGER_EDGEBLOCKSIZE; ++k) {
	  int64_t myNeighbor = tmp->edges[k].neighbor;
	  if (to == myNeighbor && k < endk) {
	    update_edge_data (G, tmp, k, to, tmp->edges[k].weight + weight, timestamp);
	    return 0;
	  }

	  if (myNeighbor < 0 || k >= endk) {
	    int64_t timefirst = readfe ( &(tmp->edges[k].timeFirst) );
	    int64_t thisEdge = tmp->edges[k].neighbor;
	    endk = tmp->high;

	    if (thisEdge < 0 || k >= endk) {
	      update_edge_data (G, tmp, k, to, weight, timestamp);
	      return 1;
	    } else if (to == thisEdge) {
	      update_edge_data (G, tmp, k, to, tmp->edges[k].weight + weight, timestamp);
	      writexf ( &(tmp->edges[k].timeFirst), timefirst);
	      return 0;
	    } else {
	      writexf ( &(tmp->edges[k].timeFirst), timefirst);
	    }
	  }
	}
	block_ptr = &(tmp->next);
      }
    }

    /* 3: Needs a new block to be inserted at end of list. */
    eb_index_t old_eb = readfe ((uint64_t *)block_ptr );
    if (!old_eb) {
      eb_index_t newBlock = new_eb (type, from);
      if (newBlock == 0) {
	writeef ((uint64_t *)block_ptr, (uint64_t)old_eb);
	return -1;
      } else {
	update_edge_data (G, ebpool_priv + newBlock, 0, to, weight, timestamp);
	ebpool_priv[newBlock].next = 0;
	push_ebs (G, 1, &newBlock);
      }
      writeef ((uint64_t *)block_ptr, (uint64_t)newBlock);
      return 1;
    }
    writeef ((uint64_t *)block_ptr, (uint64_t)old_eb);
  }
}

/** @brief Insert an undirected edge.
 *
 *  Inserts a typed, undirected edge.  First timestamp is set, if the edge is
 *  new.  Recent timestamp is updated.  Weight is set to specified value regardless.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param weight Edge weight
 *  @param timestamp Edge timestamp
 *  @return Number of edges inserted successfully
 */
MTA ("mta inline")
int
stinger_insert_edge_pair (struct stinger *G,
                          int64_t type, int64_t from, int64_t to,
                          int64_t weight, int64_t timestamp)
{
  STINGERASSERTS();

  int rtn = stinger_insert_edge (G, type, from, to, weight, timestamp);
  int rtn2 = stinger_insert_edge (G, type, to, from, weight, timestamp);

  /* Check if either returned -1 */
  if(rtn < 0 || rtn2 < 0)
    return -1;
  else
    return rtn + rtn2;
}

/** @brief Increments an undirected edge.
 *
 *  Increments the weight of a typed, undirected edge.
 *  Recent timestamp is updated.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param weight Edge weight
 *  @param timestamp Edge timestamp
 *  @return 1
 */
MTA ("mta inline")
int
stinger_incr_edge_pair (struct stinger *G,
                        int64_t type, int64_t from, int64_t to,
                        int64_t weight, int64_t timestamp)
{
  STINGERASSERTS();

  int rtn = stinger_incr_edge (G, type, from, to, weight, timestamp);
  int rtn2 = stinger_incr_edge (G, type, to, from, weight, timestamp);

  /* Check if either returned -1 */
  if(rtn < 0 || rtn2 < 0)
    return -1;
  else
    return rtn + rtn2;
}

/** @brief Removes a directed edge.
 *
 *  Remove a typed, directed edge.
 *  Note: Do not call this function concurrently with the same source vertex,
 *  even for different edge types.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @return 1 on success, 0 if the edge is not found.
 */
MTA ("mta inline") MTA ("mta serial")
int
stinger_remove_edge (struct stinger *G,
                     int64_t type, int64_t from, int64_t to)
{
  if(from == to)
    return -1;

  /* Do *NOT* call this concurrently with different edge types. */
  STINGERASSERTS ();

  struct curs curs;
  struct stinger_eb *tmp;
  struct stinger_eb *ebpool_priv = ebpool;

  curs = etype_begin (&G->LVA[from], type);

  for (tmp = ebpool_priv + curs.eb; tmp != ebpool_priv; tmp = ebpool_priv + readff((uint64_t *)&tmp->next)) {
    if(type == tmp->etype) {
      size_t k, endk;
      endk = tmp->high;

      for (k = 0; k < endk; ++k) {
	if (to == tmp->edges[k].neighbor) {
	  int64_t weight = readfe (&(tmp->edges[k].weight));
	  if(to == tmp->edges[k].neighbor) {
	    update_edge_data (G, tmp, k, ~to, weight, 0);
	    return 1;
	  } else {
	    writeef((uint64_t *)&(tmp->edges[k].weight), (uint64_t)weight);
	  }
	  return 0;
	}
      }
    }
  }
}

/** @brief Removes an undirected edge.
 *
 *  Remove a typed, undirected edge.
 *  Note: Do not call this function concurrently with the same source vertex,
 *  even for different edge types.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @return 1 on success, 0 if the edge is not found.
 */
MTA ("mta inline")
int
stinger_remove_edge_pair (struct stinger *G,
                          int64_t type, int64_t from, int64_t to)
{
  STINGERASSERTS();

  int rtn = stinger_remove_edge (G, type, from, to);
  int rtn2 = stinger_remove_edge (G, type, to, from);

  /* Check if either returned -1 */
  if(rtn < 0 || rtn2 < 0)
    return -1;
    else
    return rtn + rtn2;
}

MTA("mta parallel default")

/** @brief Initializes an empty STINGER with a graph in CSR format.
 *
 *  Takes an edge list in CSR format, with weight and timestamps, and initializes
 *  an empty STINGER based on the input graph.  All edges being ingested must be
 *  of a single edge type.
 *
 *  @param G STINGER data structure
 *  @param nv Number of vertices
 *  @param etype Edge type
 *  @param off_in Array of length nv containing the adjacency offset for each vertex
 *  @param phys_adj_in Array containing the destination vertex of each edge
 *  @param weight_in Array containing integer weight of each edge
 *  @param ts_in Array containing recent timestamp of edge edge (or NULL)
 *  @param first_ts_in Array containing the first timestamp of each edge (or NULL)
 *  @param single_ts Value for timestamps if either of the above is NULL
 *  @return Void
 */
MTA ("mta inline")
void
stinger_set_initial_edges (struct stinger *G,
                           const size_t nv,
                           const int64_t etype,
                           const int64_t * off_in,
                           const int64_t * phys_adj_in,
                           const int64_t * weight_in,
                           const int64_t * ts_in,
                           const int64_t * first_ts_in,
                           const int64_t single_ts
                           /* if !ts or !first_ts */ )
{
  const int64_t *restrict off = off_in;
  const int64_t *restrict phys_adj = phys_adj_in;
  const int64_t *restrict weight = weight_in;
  const int64_t *restrict ts = ts_in;
  const int64_t *restrict first_ts = first_ts_in;
  struct stinger_vb *restrict LVA;

  size_t nblk_total = 0;
  size_t *restrict blkoff;
  eb_index_t *restrict block;

  assert (G);
  LVA = G->LVA;

  blkoff = xcalloc (nv + 1, sizeof (*blkoff));
  OMP ("omp parallel for")
    for (int64_t v = 0; v < nv; ++v) {
      const int64_t deg = off[v + 1] - off[v];
      blkoff[v + 1] = (deg + STINGER_EDGEBLOCKSIZE - 1) / STINGER_EDGEBLOCKSIZE;
    }

  for (int64_t v = 2; v <= nv; ++v)
    blkoff[v] += blkoff[v - 1];
  nblk_total = blkoff[nv];

  block = xcalloc (nblk_total, sizeof (*block));
  OMP ("omp parallel for")
    MTA ("mta assert nodep") MTASTREAMS ()
    for (int64_t v = 0; v < nv; ++v) {
      const int64_t from = v;
      const int64_t deg = off[v + 1] - off[v];
      if (deg)
        stinger_int64_fetch_add (&LVA[from].outDegree, deg);
    }

  new_blk_ebs (&block[0], G, nv, blkoff, etype);

  /* XXX: AUGH! I cannot find what really is blocking parallelization. */
  MTA ("mta assert parallel")
    MTA ("mta dynamic schedule")
    OMP ("omp parallel for")
    MTA ("mta assert noalias *block")
    MTA ("mta assert nodep *block")
    MTA ("mta assert noalias *G")
    MTA ("mta assert nodep *G")
    MTA ("mta assert nodep *phys_adj")
    MTA ("mta assert nodep *LVA")
    MTA ("mta assert nodep *blkoff")
    for (int64_t v = 0; v < nv; ++v) {
      const size_t nextoff = off[v + 1];
      size_t kgraph = off[v];
      MTA ("mta assert may reorder kgraph") int64_t from;

      from = v;

      // Need to assert this loop is not to be parallelized!
      MTA ("mta assert nodep *block")
        MTA ("mta assert nodep *LVA")
        for (size_t kblk = blkoff[v]; kblk < blkoff[v + 1]; ++kblk) {
          size_t n_to_copy, voff;
          struct stinger_edge *restrict edge;
          struct stinger_eb *restrict eb;
          int64_t tslb = INT64_MAX, tsub = 0;

          //voff = stinger_size_fetch_add (&kgraph, STINGER_EDGEBLOCKSIZE);
          {
            voff = kgraph;
            kgraph += STINGER_EDGEBLOCKSIZE;
          }

          n_to_copy = STINGER_EDGEBLOCKSIZE;
          if (voff + n_to_copy >= nextoff)
            n_to_copy = nextoff - voff;

          eb = ebpool + block[kblk];
          edge = &eb->edges[0];

          /* XXX: remove the next two asserts once the outer is unlocked. */
          MTA ("mta assert nodep")
            MTASTREAMS ()MTA ("mta assert nodep *phys_adj")
            MTA ("mta assert nodep *edge")
            MTA ("mta assert nodep *LVA")
            for (size_t i = 0; i < n_to_copy; ++i) {
              const int64_t to = phys_adj[voff + i];
#if defined(__MTA__)
              int_fetch_add (&LVA[to].inDegree, 1);
#else
              // Ugh. The MTA compiler can't cope with the inlining.
              stinger_int64_fetch_add (&LVA[to].inDegree, 1);
#endif
              /* XXX: The next statements block parallelization
                 of the outer loop. */
              edge[i].neighbor = to;
              edge[i].weight = weight[voff + i];
              edge[i].timeRecent = ts ? ts[voff + i] : single_ts;
              edge[i].timeFirst = first_ts ? first_ts[voff + i] : single_ts;
              //assert (edge[i].timeRecent >= edge[i].timeFirst);
            }

          if (ts || first_ts) {
            /* XXX: remove the next two asserts once the outer is unlocked. */
            MTA ("mta assert nodep")
              MTASTREAMS ()MTA ("mta assert nodep *edge")
              for (size_t i = 0; i < n_to_copy; ++i) {
                if (edge[i].timeFirst < tslb)
                  tslb = edge[i].timeFirst;
                if (edge[i].timeRecent > tsub)
                  tsub = edge[i].timeRecent;
              }
          } else
            tslb = tsub = single_ts;

          eb->smallStamp = tslb;
          eb->largeStamp = tsub;
          eb->numEdges = n_to_copy;
          eb->high = n_to_copy;
        }


      /* At this point, block[blkoff[v]] is the head of a linked
         list holding all the blocks of edges of EType from vertex
         v.  Insert into the graph.  */

      if (blkoff[v] != blkoff[v + 1])
        LVA[from].edges = block[blkoff[v]];
      else
        LVA[from].edges = 0;
    }

  /* Insert into the edge type array */
  push_ebs (G, nblk_total, block);

  free (block);
  free (blkoff);
}

/** @brief Copy typed incoming adjacencies of a vertex into a buffer
 *
 *  For a given edge type, adjacencies of the specified vertex are copied into
 *  the user-provided buffer up to the length of the buffer.  These are the
 *  incoming edges for which a vertex is a destination.  Note that this operation
 *  may be very expensive on most platforms.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param v Source vertex ID
 *  @param outlen Number of adjacencies copied
 *  @param out Buffer to hold adjacencies
 *  @param max_outlen Length of out[] and recent[]
 *  @return Void
 */
MTA ("mta inline")
void
stinger_gather_typed_predecessors (const struct stinger *G,
				   int64_t type,
				   int64_t v,
				   size_t * outlen,
                                   int64_t * out,
                                   size_t max_outlen)
{
  size_t kout = 0;

  assert (G);

  STINGER_PARALLEL_FORALL_EDGES_BEGIN(G, type) {
    const int64_t u = STINGER_EDGE_SOURCE;
    const int64_t v = STINGER_EDGE_DEST;
    if (v >= 0) {
      size_t where = stinger_size_fetch_add (&kout, 1);
      if (where < max_outlen) {
        out[where] = u;
      }
    }
  } STINGER_PARALLEL_FORALL_EDGES_END();

  *outlen = kout;               /* May be longer than max_outlen. */
}

/** @brief Copy adjacencies of a vertex into a buffer with optional metadata
 *
 *  Adjacencies of the specified vertex are copied into the user-provided buffer(s) 
 *  up to the length of the buffer(s) specified by max_outlen.  All buffers should 
 *  be at least max_outlen or NULL.
 *
 *  @param G The STINGER data structure
 *  @param v Source vertex ID
 *  @param outlen Number of adjacencies copied
 *  @param out Buffer to hold adjacencies
 *  @param weight OPTIONAL Buffer to hold edge weights
 *  @param timefirst OPTIONAL Buffer to hold first timestamps
 *  @param timerecent OPTIONAL Buffer to hold recent timestamps
 *  @param type OPTIONAL Buffer to hold edge types
 *  @param max_outlen Length of out[] and any optional buffers provided
 *  @return Void
 */
MTA ("mta inline")
void
stinger_gather_successors (const struct stinger *G,
                                          int64_t v,
                                          size_t * outlen,
                                          int64_t * out,
					  int64_t * weight,
					  int64_t * timefirst,
                                          int64_t * timerecent,
					  int64_t * type,
                                          size_t max_outlen)
{
  size_t kout = 0;

  assert (G);

  STINGER_PARALLEL_FORALL_EDGES_OF_VTX_BEGIN(G, v) {
    const int64_t n = STINGER_EDGE_DEST;
    const int64_t w = STINGER_EDGE_WEIGHT;
    const int64_t tf = STINGER_EDGE_TIME_FIRST;
    const int64_t tr = STINGER_EDGE_TIME_RECENT;
    const int64_t t = STINGER_EDGE_TYPE;
    if (n >= 0) {
      size_t where = stinger_size_fetch_add (&kout, 1);
      if (where < max_outlen) {
        out[where] = n;
        if(weight) weight[where] = w;
        if(timefirst) timefirst[where] = tf;
        if(timerecent) timerecent[where] = tr;
        if(type) type[where] = t;
      }
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_VTX_END();

  *outlen = kout;               /* May be longer than max_outlen. */
}

/** @brief Copy typed adjacencies of a vertex into a buffer
 *
 *  For a given edge type, adjacencies of the specified vertex are copied into
 *  the user-provided buffer up to the length of the buffer.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param v Source vertex ID
 *  @param outlen Number of adjacencies copied
 *  @param out Buffer to hold adjacencies
 *  @param max_outlen Length of out[] and recent[]
 *  @return Void
 */
MTA ("mta inline")
void
stinger_gather_typed_successors (const struct stinger *G, int64_t type,
                                 int64_t v, size_t * outlen,
                                 int64_t * out, size_t max_outlen)
{
  size_t kout = 0;

  assert (G);

  if (v < 0) {
    *outlen = 0;
    return;
  }

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G, type, v) {
    size_t where = stinger_size_fetch_add (&kout, 1);
    if (where < max_outlen)
      out[where] = STINGER_EDGE_DEST;
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();

  *outlen = kout;               /* May be longer than max_outlen. */
}

/** @brief Determines if a vertex has an edge of a given type
 *
 *  Searches the adjacencies of a specified vertex for an edge of the given type.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @return 1 if found; 0 otherwise
 */
MTA ("mta inline")
int
stinger_has_typed_successor (const struct stinger *G,
    int64_t type, int64_t from, int64_t to)
{
  STINGERASSERTS ();

  int rtn = 0;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      stinger_int_fetch_add(&rtn, 1);
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return (rtn > 0 ? 1 : 0); 
}

/** @brief Get the weight of a given edge.
 *
 *  Finds a specified edge of a given type by source and destination vertex ID
 *  and returns the current edge weight.  Remember, edges may have different
 *  weights in different directions.
 *
 *  @param G The STINGER data structure
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param type Edge type
 *  @return Edge weight
 */
MTA ("mta inline")
int64_t
stinger_edgeweight (const struct stinger * G,
                    int64_t from, int64_t to, int64_t type)
{
  STINGERASSERTS ();

  int rtn = 0;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      rtn = STINGER_EDGE_WEIGHT;
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return rtn;
}

/* TODO revisit this function
 * XXX how to handle in parallel?
 * BUG block meta not handled
 * */
/** @brief Set the weight of a given edge.
 *
 *  Finds a specified edge of a given type by source and destination vertex ID
 *  and sets the current edge weight.  Remember, edges may have different
 *  weights in different directions.
 *
 *  @param G The STINGER data structure
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param type Edge type
 *  @param weight Edge weight to set
 *  @return 1 on success; 0 otherwise
 */
MTA ("mta inline")
int
stinger_set_edgeweight (struct stinger *G,
                        int64_t from, int64_t to, int64_t type,
                        int64_t weight)
{
  STINGERASSERTS ();

  int rtn = 0;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      STINGER_EDGE_WEIGHT = weight;
      rtn = 1;
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return rtn;
}

/** @brief Get the first timestamp of a given edge.
 *
 *  Finds a specified edge of a given type by source and destination vertex ID
 *  and returns the first timestamp.  Remember, edges may have different
 *  timestamps in different directions.
 *
 *  @param G The STINGER data structure
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param type Edge type
 *  @return First timestamp of edge; -1 if edge does not exist
 */
MTA ("mta inline")
int64_t
stinger_edge_timestamp_first (const struct stinger * G,
                              int64_t from, int64_t to, int64_t type)
{
  STINGERASSERTS ();

  int rtn = -1;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      rtn = STINGER_EDGE_TIME_FIRST;
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return rtn;
}

/** @brief Get the recent timestamp of a given edge.
 *
 *  Finds a specified edge of a given type by source and destination vertex ID
 *  and returns the recent timestamp.  Remember, edges may have different
 *  timestamps in different directions.
 *
 *  @param G The STINGER data structure
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param type Edge type
 *  @return Recent timestamp of edge; -1 if edge does not exist
 */
MTA ("mta inline")
int64_t
stinger_edge_timestamp_recent (const struct stinger * G,
                               int64_t from, int64_t to, int64_t type)
{
  STINGERASSERTS ();

  int rtn = -1;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      rtn = STINGER_EDGE_TIME_RECENT;
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return rtn;
}

/* TODO revisit this function
 * XXX how to handle in parallel?
 * BUG block meta not handled
 * */
/** @brief Update the recent timestamp of an edge
 *
 *  Finds a specified edge of a given type by source and destination vertex ID
 *  and updates the recent timestamp to the specified value.
 *
 *  @param G The STINGER data structure
 *  @param from Source vertex ID
 *  @param to Destination vertex ID
 *  @param type Edge type
 *  @param timestamp Timestamp to set recent
 *  @return 1 on success; 0 if failure
 */
MTA ("mta inline")
int
stinger_edge_touch (struct stinger *G,
    int64_t from, int64_t to, int64_t type,
    int64_t timestamp)
{
  STINGERASSERTS ();

  int rtn = 0;

  STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G,type,from) {
    if (STINGER_EDGE_DEST == to) {
      STINGER_EDGE_TIME_RECENT = timestamp;
      rtn = 1;
    }
  } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();
  return rtn;
}

/* TODO - add doxygen and insert into stinger.h */
MTA ("mta inline") int64_t
stinger_count_outdeg (struct stinger * G, int64_t v)
{
  int64_t nactive_edges = 0, neb = 0;

  STINGER_FORALL_EB_BEGIN (G, v, eb) {
    const size_t eblen = stinger_eb_high (eb);
    OMP ("omp parallel for")
      MTA ("mta assert nodep")
      MTASTREAMS ()
      for (size_t ek = 0; ek < eblen; ++ek) {
        if (!stinger_eb_is_blank (eb, ek))
          stinger_int64_fetch_add (&nactive_edges, 1);
      }
  }
  STINGER_FORALL_EB_END ();

  return nactive_edges;
}

/**
* @brief Sorts a batch of edge insertions and deletions.
* 
* Takes an array of edge insertions and deletions and sorts them.  The array is
* packed with <source, destination> pairs, such that actions[2*i] = source vertex
* ID and actions[2*i+1] = destination vertex ID.  Bit-wise complement the source
* and destination vertex IDs to indicate a deletion.  This function will create
* an undirected actions list (i.e. create the reverse edge action).
*
* The output is similarly packed, sorted by source vertex ID first, then by
* destination vertex ID such that deletions are contiguous before insertions.
* For vertex i, deletions incident on i start at deloff[i] and end at insoff[i].
* Insertions incident on i start at insoff[i] and end at deloff[i+1].
*
* @param nactions Number of edge actions
* @param actions The packed array of edge insertions and deletions
* @param insoff For each incident vertex in the batch, the offset into act[] of
*   the first edge insertion
* @param deloff For each incident vertex in the batch, the offset into act[] of
*   the first edge deletion
* @param act Sorted array of edge insertions and deletions
*
* @return The number of incident vertices in the batch (the size of insoff[] and deloff[])
*/
int64_t
stinger_sort_actions (int64_t nactions, int64_t * actions,
                      int64_t * insoff, int64_t * deloff, int64_t * act)
{
  int64_t head = 0;
  int64_t actlen;
  int64_t n;
  int64_t *actk = xmalloc ((2 * 2 * nactions + 1) * sizeof (*actk));

  OMP("omp parallel") {
  /* Copy & make i positive if necessary.
     Negative j still implies deletion. */
  MTA ("mta assert nodep")
    MTA ("mta block schedule")
    OMP ("omp for")
    for (int64_t k = 0; k < nactions; k++) {
      const int64_t i = actions[2 * k];
      const int64_t j = actions[2 * k + 1];

      if (i != j) {
        int64_t where = stinger_int64_fetch_add (&head, 2);
        assert (where < 2 * nactions);

        act[2 * where] = (i < 0 ? -i - 1 : i);
        act[2 * where + 1] = j;
        act[2 * (where + 1)] = (j < 0 ? -j - 1 : j);
        act[2 * (where + 1) + 1] = i;
      }
    }

  OMP("omp single") {
  actlen = head;
#if !defined(__MTA__)
  qsort (act, actlen, 2 * sizeof (act[0]), i2cmp); 
  //radix_sort_pairs (act, actlen<<1, 6);
#else
  bucket_sort_pairs (act, actlen);
#endif

  actk[0] = 0;
  n = 1;
  }

  /* Find local indices... */
  MTA ("mta assert nodep")
    OMP ("omp for")
    for (int64_t k = 1; k < actlen; k++) {
      if (act[2 * k] != act[2 * (k - 1)]) {
        stinger_int64_fetch_add (&n, 1);
        actk[k] = 1;
      } else
        actk[k] = 0;
    }

  prefix_sum (actlen, actk);
  //XXX assert(actk[actlen-1] == n-1);

  MTA ("mta assert nodep")
    OMP ("omp for")
    for (int64_t k = 0; k <= n; k++)
      deloff[k] = 0;

  OMP ("omp for")
    MTA ("mta assert nodep")
    for (int64_t k = 0; k < actlen; k++)
      stinger_int64_fetch_add (&deloff[actk[k] + 1], 1);

  prefix_sum (n+1, deloff);
  assert (deloff[n] == actlen);

  MTA ("mta assert nodep")
    MTA ("mta loop norestructure")
    OMP ("omp for")
    for (int64_t k = 0; k < n; k++) {
      int off;
      const int endoff = deloff[k + 1];

      for (off = deloff[k]; off < endoff; off++)
        if (act[2 * off + 1] >= 0)
          break;
      insoff[k] = off;
      assert (insoff[k] == deloff[k + 1]
              || act[2 * insoff[k]] == act[2 * deloff[k]]);
      for (off = deloff[k]; off < insoff[k]; off++) {
        assert (act[2 * off + 1] < 0);
        act[2 * off + 1] = -act[2 * off + 1] - 1;
      }
    }
  }

  free (actk);

  return n;                     /* number of unique vertices, elements in insoff[], deloff[] */
}

/** @brief Removes all the edges in the graph of a given type.
 *
 *  Traverses all edge blocks of a particular type in parallel and erases all
 *  edges of the specified type.  Blocks are available immediately for reuse.
 *
 *  @param G The STINGER data structure
 *  @param type The edge type of edges to delete
 *  @return Void
 */
void
stinger_remove_all_edges_of_type (struct stinger *G, int64_t type)
{
  /* TODO fix bugs here */
  MTA("mta assert parallel")
  OMP("omp parallel for")
  for (uint64_t p = 0; p < G->ETA[type].high; p++) {
    struct stinger_eb *current_eb = ebpool + G->ETA[type].blocks[p];
    int64_t thisVertex = current_eb->vertexID;
    int64_t high = current_eb->high;
    struct stinger_edge * edges = current_eb->edges;

    int64_t removed = 0;
    for (uint64_t i = 0; i < high; i++) {
      int64_t neighbor = edges[i].neighbor;
	if (neighbor >= 0) {
	removed++;
	assert(neighbor >= 0);
	stinger_int64_fetch_add (&(G->LVA[neighbor].inDegree), -1);
	edges[i].neighbor = -1;
      }
    }
    stinger_int64_fetch_add (&(G->LVA[thisVertex].outDegree), -removed);
    current_eb->high = 0;
    current_eb->numEdges = 0;
    current_eb->smallStamp = INT64_MAX;
    current_eb->largeStamp = INT64_MIN;
  }
}

/** @brief Checkpoint a STINGER data structure to disk.
 *
 *  Creates a directory if needed and stores the graph as individual adjacency files in directory
 *  for parallelism / no need for extra buffers.
 *  Format:
 *  endian_check -> 64-bit endianness file
 *  nv -> 64-bit max vertex ID
 *  0 .. maxVtx -> 64-bit outdegree followed by adjacencies of vertex 0 stored as 5-tuples 
 *     (type, dest, weight, timefirst, timerecent)
 *
 *  @param S The STINGER data structure
 *  @param maxVtx Largest logical vertex ID + 1
 *  @param stingerfiledir String containing the name of the directory to store the checkpoint files in
 *  @return 0 on success; -1 if error
 */
int
stinger_save_to_file (struct stinger * S, uint64_t maxVtx, const char * stingerfiledir)
{
   /* TODO: 
   * - XMT testing / code
   * - Performance / size tweaks (this is a rough first run)
   */

#if !defined(__MTA__)
  if(stinger_consistency_check(S, maxVtx)) {
    fprintf(stderr, "WARNING: Stinger is inconsistent before writing to disk.\n");
  }

  mkdir(stingerfiledir, S_IRWXU | S_IRWXG | S_IROTH);

  if(chdir(stingerfiledir)) {
    fprintf(stderr, "ERROR: Directory creating / cd Failed.\n");
    return -1;
  }

  struct dirent **names;
  int64_t n = scandir(".", &names, 0, alphasort);
  if(n >= 0) {
    while(n--) {
      if(names[n]->d_type == DT_REG)
	unlink(names[n]->d_name);
      free(names[n]);
    }
  }
  free(names);

  {
    const int64_t endian_check = 0x1234ABCD;
    FILE * f_endianness = fopen("endian_check", "w"); 
    if(!f_endianness) {
      fprintf(stderr, "ERROR: Creating endianness failed.\n");
      return -1;
    }
    if(1 != fwrite(&endian_check, sizeof(int64_t), 1, f_endianness)) {
      fprintf(stderr, "ERROR: Writing endianness failed.\n");
      return -1;
    }
    fclose(f_endianness);
  }

  {
    FILE * f_nv = fopen("nv", "w"); 
    if(!f_nv || 1 != fwrite(&maxVtx, sizeof(int64_t), 1, f_nv)) {
      fprintf(stderr, "ERROR: Writing NV failed.\n");
      return -1;
    }
    fclose(f_nv);
  }

  uint64_t error = 0;

  OMP("omp parallel for")
  for(uint64_t i = 0; i < maxVtx; i++) {
    uint64_t outDegree = stinger_outdegree(S, i);
    if(outDegree) {
      char filename[256];
      sprintf(filename, "%lu", i);
      FILE * f_curvtx = fopen(filename, "w");
      if(!f_curvtx || 1 != fwrite(&outDegree, sizeof(uint64_t), 1, f_curvtx)) {
	fprintf(stderr, "ERROR: Writing %lu outdegree failed.\n", i);
	error = -1;
      } else {
	uint64_t value = stinger_vweight(S, i);
	if(1 != fwrite(&value, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Writing %lu vweight failed.\n", i);
	  error = -1;
	}
	value = stinger_vtype(S, i);
	if(1 != fwrite(&value, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Writing %lu vtype failed.\n", i);
	  error = -1;
	}
	STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, i) {
	  if(1 != fwrite(&(STINGER_EDGE_TYPE), sizeof(uint64_t), 1, f_curvtx)) {
	    fprintf(stderr, "ERROR: Writing %lu type failed.\n", i);
	    error = -1;
	  }
	  if(1 != fwrite(&(STINGER_EDGE_DEST), sizeof(uint64_t), 1, f_curvtx)) {
	    fprintf(stderr, "ERROR: Writing %lu dest failed.\n", i);
	    error = -1;
	  }
	  if(1 != fwrite(&(STINGER_EDGE_WEIGHT), sizeof(uint64_t), 1, f_curvtx)) {
	    fprintf(stderr, "ERROR: Writing %lu weight failed.\n", i);
	    error = -1;
	  }
	  if(1 != fwrite(&(STINGER_EDGE_TIME_FIRST), sizeof(uint64_t), 1, f_curvtx)) {
	    fprintf(stderr, "ERROR: Writing %lu timefirst failed.\n", i);
	    error = -1;
	  }
	  if(1 != fwrite(&(STINGER_EDGE_TIME_RECENT), sizeof(uint64_t), 1, f_curvtx)) {
	    fprintf(stderr, "ERROR: Writing %lu timerecent failed.\n", i);
	    error = -1;
	  }
	} STINGER_FORALL_EDGES_OF_VTX_END();
	fclose(f_curvtx);
      }
    }
  }
  if(chdir("..")) {
    fprintf(stderr, "ERROR: Returning to parent.\n");
    error = -1;
  }

  return error;
#endif
}

#define EBSTACK_NEW(X, V) \
  int64_t ebstack_top = 0; \
  int64_t ebstack_size = X;\
  eb_index_t * ebstack = xmalloc(X * sizeof(eb_index_t)); \
  new_ebs(ebstack, X, 0, V);

#define EBSTACK_POP(Y, V) \
  if(ebstack_top == ebstack_size) { \
    push_ebs(*S, ebstack_size, ebstack); \
    ebstack_size = 1; \
    ebstack_top = 0;  \
    new_ebs(ebstack, 1, 0, V);  \
  } \
  Y = ebstack[ebstack_top++];

#define EBSTACK_FREE() \
  push_ebs(*S, ebstack_size, ebstack); \
  free(ebstack);

/* TODO XXX BUG push_ebs won't work with the ebstack correctly here since
 * ebs in the stack can be different types */
/** @brief Restores a STINGER checkpoint from disk.
 *
 *  Reads STINGER checkpoint files from disk and reconstructs the STINGER data structure
 *  according to their contents.
 *
 *  @param stingerfiledir Directory name containing the STINGER checkpoint files
 *  @param S STINGER data structure to restore into
 *  @param maxVtx Largest logical vertex ID + 1
 *  @return 0 on success; -1 if error 
 */
int
stinger_open_from_file (const char * stingerfiledir, struct stinger ** S, uint64_t * maxVtx)
{

#if !defined(__MTA__)
  *S = NULL;
  if(chdir(stingerfiledir)) {
    fprintf(stderr, "ERROR: Changing to directory.\n");
    *S = NULL;
    return -1;
  }
  int64_t endian_check = 0x1234ABCD;
  {
    FILE * f_endianness = fopen("endian_check", "r");
    int64_t endian_input;
    if(!f_endianness || 1 != fread(&endian_input, sizeof(int64_t), 1, f_endianness)) {
      fprintf(stderr, "ERROR: Can't read endianness file.\n");
      *S = NULL;
      return -1;
    }
    endian_check = endian_check != endian_input;
  }

  {
    FILE * f_nv = fopen("nv", "r");
    if(!f_nv || 1 != fread(maxVtx, sizeof(int64_t), 1, f_nv)) {
      fprintf(stderr, "ERROR: Can't read maxvtx file.\n");
      *S = NULL;
      return -1;
    }
    if(endian_check) {
      *maxVtx = bs64(*maxVtx);
    }
  }

  *S = stinger_new();

  int error = 0;
  OMP("omp parallel for")
  for(uint64_t i = 0; i < *maxVtx; i++) {
    char filename[256];
    sprintf(filename, "%lu", i);
    FILE * f_curvtx = fopen(filename, "r");
    if(f_curvtx) {
      uint64_t outDegree;
      if(1 != fread(&outDegree, sizeof(int64_t), 1, f_curvtx)) {
	fprintf(stderr, "ERROR: Can't read %lu file.\n", i);
	error = -1;
      }
      if(endian_check) {
	outDegree = bs64(outDegree);
      }
      uint64_t value = 0;
      if(1 != fread(&value, sizeof(int64_t), 1, f_curvtx)) {
	fprintf(stderr, "ERROR: Can't read %lu vweight file.\n", i);
	error = -1;
      }
      if(endian_check) {
	value = bs64(value);
      }
      stinger_set_vweight(*S, i, value);
      if(1 != fread(&value, sizeof(int64_t), 1, f_curvtx)) {
	fprintf(stderr, "ERROR: Can't read %lu vtype file.\n", i);
	error = -1;
      }
      if(endian_check) {
	value = bs64(value);
      }
      stinger_set_vtype(*S, i, value);
      EBSTACK_NEW((outDegree + STINGER_EDGEBLOCKSIZE - 1) / STINGER_EDGEBLOCKSIZE, i);
      int64_t type, dest, weight, tfirst, trecent;
      int read;
      EBSTACK_POP((*S)->LVA[i].edges, i);
      struct stinger_eb * curblock = ebpool + (*S)->LVA[i].edges;
      while((read = fread(&type, sizeof(uint64_t), 1, f_curvtx)) != 0) {
	if(curblock->high == STINGER_EDGEBLOCKSIZE) {
	  EBSTACK_POP(eb_index_t nextblock, i);
	  curblock->next = nextblock;
	  curblock = ebpool + nextblock;
	}

	if(endian_check) {
	  type = bs64(type);
	}
	// handle type
	if(!curblock->numEdges) {
	  curblock->etype = type;
	} else {
	  if(curblock->etype != type) {
	    EBSTACK_POP(eb_index_t nextblock, i);
	    curblock->next = nextblock;
	    curblock = ebpool + nextblock;
	    curblock->etype = type;
	  }
	}

	if(1 != fread(&dest, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Can't read %lu. Unexpected end\n", i);
	  error = -1;
	}
	if(endian_check) {
	  dest = bs64(dest);
	  curblock->etype = type;
	} else {
	  if(curblock->etype != type) {
	    EBSTACK_POP(eb_index_t nextblock, i);
	    curblock->next = nextblock;
	    curblock = ebpool + nextblock;
	    curblock->etype = type;
	  }
	}

	if(1 != fread(&dest, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Can't read %lu. Unexpected end\n", i);
	  error = -1;
	}
	if(endian_check) {
	  dest = bs64(dest);
	}
	if(1 != fread(&weight, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Can't read %lu. Unexpected end\n", i);
	  error = -1;
	}
	if(endian_check) {
	  weight = bs64(weight);
	}
	if(1 != fread(&tfirst, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Can't read %lu. Unexpected end\n", i);
	  error = -1;
	}
	if(endian_check) {
	  tfirst= bs64(tfirst);
	}
	if(1 != fread(&trecent, sizeof(uint64_t), 1, f_curvtx)) {
	  fprintf(stderr, "ERROR: Can't read %lu. Unexpected end\n", i);
	  error = -1;
	}
	if(endian_check) {
	  trecent = bs64(trecent);
	}
	stinger_int64_fetch_add(&((*S)->LVA[i].outDegree), 1);
	stinger_int64_fetch_add(&((*S)->LVA[dest].inDegree), 1);
	curblock->edges[curblock->high].neighbor = dest;
	curblock->edges[curblock->high].weight = weight;
	curblock->edges[curblock->high].timeFirst = tfirst;
	curblock->edges[curblock->high].timeRecent = trecent;
	if(curblock->smallStamp > tfirst)
	  curblock->smallStamp = tfirst;
	if(curblock->largeStamp < trecent)
	  curblock->largeStamp = trecent;
	curblock->numEdges++;
	curblock->high++;

      }
      fclose(f_curvtx);
      EBSTACK_FREE();
    }
  }

  if(chdir("..")) {
    fprintf(stderr, "ERROR: Returning to parent.\n");
    error = -1;
  }
  return error;
  if(error)
    return error;
  else
    return stinger_consistency_check(*S, *maxVtx);

#endif
}


#include "stinger-utils.h"
#include "xmalloc.h"
#include "stinger-atomics.h"
#include "stinger.h"

MTA ("mta parallel off") MTA ("mta expect parallel context")
/** @brief DEPRECATED Remove and insert edges incident on a common source vertex.
 *
 *  DEPRECATED
 *  For a given source vertex and edge type, removes a list of edges given in
 *  an array.  Then inserts a list of edges given in a second array with weight
 *  and a timestamp.
 *  Note: Do not call this function concurrently with the same source vertex,
 *  even for different edge types.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param from Source vertex ID
 *  @param nremove Number of edges to remove
 *  @param remove Edge adjacencies to remove
 *  @param ninsert Number of edges to insert
 *  @param insert Edge adjacencies to insert
 *  @param weight Edge adjacency weights to insert
 *  @param timestamp Timestamp for all edge insertions
 *  @return 0 on success; number of failures otherwise
 */
int
stinger_remove_and_insert_edges (struct stinger *G,
                                 int64_t type, int64_t from,
                                 int64_t nremove, int64_t * remove,
                                 int64_t ninsert, int64_t * insert,
                                 int64_t * weight, int64_t timestamp)
{
  /* Do *NOT* call this concurrently with the same from value, even
     with different types. */
  /*
    First pass:
    - Remove vertices.
    - Scan for insertion slots.
    - Scan for existing vertices being inserted.
    Second:
    - Insert into saved locations.
  */

  int64_t nrem = 0;
  struct curs curs;
  struct stinger_eb *tmp;
  eb_index_t *prev_loc;
  struct stinger_eb **has_slot = NULL;
  int *kslot = NULL;
  int64_t nslot = 0, ninsert_remaining = ninsert;

#if defined(__MTA__)
  if (ninsert < 4) {
    for (int64_t k = 0; k < 4 && ninsert > 0; k++) {
      if (insert[k] >= 0)
	nrem += stinger_insert_edge (G, type, from, insert[k],
				    (weight ? weight[0] : 1), timestamp);
      ninsert--;
    }
  }
#else
  if (1 == ninsert) {
    if (insert[0] >= 0)
      nrem += stinger_insert_edge (G, type, from, insert[0],
                                  (weight ? weight[0] : 1), timestamp);
    ninsert = 0;
  }
#endif

#if defined(__MTA__)
  if (nremove < 2) {
    for (int64_t k = 0; k < 2 && nremove > 0; k++) {
      if (remove[k] >= 0)
	nrem += stinger_remove_edge (G, type, from, remove[k]);
      nremove--;
    }
  }
#else
  if (1 == nremove) {
    if (remove[0] >= 0)
      nrem = stinger_remove_edge (G, type, from, remove[0]);
    nremove = 0;
  }
#endif

  /* Sort for binary search */
  if (nremove) {
    qsort (remove, nremove, sizeof (*remove), i64_cmp);
    while (nremove && *remove < 0) {
      ++remove;
      --nremove;
    }
  }
  if (ninsert) {
    has_slot = xmalloc (ninsert * sizeof (*has_slot));
    kslot = xmalloc (ninsert * sizeof (*kslot));
    qsort (insert, ninsert, sizeof (*insert), i64_cmp);
    while (ninsert && *insert < 0) {
      ++insert;
      --ninsert;
    }
  }

  if (!nremove && !ninsert) {
    return nrem;
  }

  curs = etype_begin (&G->LVA[from], type);
  prev_loc = curs.loc;

  for (tmp = ebpool + curs.eb; tmp != ebpool; tmp = ebpool + tmp->next) {
    if(tmp->etype == type) {
      size_t k, endk;
      size_t found_k = STINGER_EDGEBLOCKSIZE;
      size_t found_slot = STINGER_EDGEBLOCKSIZE;
      prev_loc = &tmp->next;
      endk = tmp->high;
      ninsert_remaining = ninsert;

	for (k = 0; k < endk; ++k) {
	  const int64_t w = tmp->edges[k].neighbor;
	  int64_t off;

	  if (w >= 0) {
	    /* Removal first.  If the vertex is removed,
	       save its slot for reuse. */
	    if (nremove) {
	      off = find_in_sorted (w, nremove, remove);
	      if (off >= 0) {
		int64_t which = stinger_int64_fetch_add (&nslot, 1);
		if (which < ninsert_remaining) {    /* Can be racy. */
		  has_slot[which] = tmp;
		  kslot[which] = k;
		}
		update_edge_data(G, tmp, k, ~w, 0, timestamp);
		stinger_int64_fetch_add (&nrem, 1);
	      }
	    }
	    if (ninsert) {
	      /* Check if this is to be inserted. */
	      off = find_in_sorted (w, ninsert, insert);
	      while (off >= 0) {
		stinger_int64_fetch_add (&ninsert_remaining, -1);
		insert[off] = ~insert[off];
		update_edge_data(G, tmp, k, w, (weight ? weight[off] : 1), timestamp);
		qsort (insert, ninsert, sizeof (*insert), i64_cmp);  /* Must maintain sorted order here. */
		off = find_in_sorted (w, ninsert, insert);   /* Gotta check if there's another one. */
	      }
	    }
	  } else if (ninsert_remaining) {
	    /* Empty slot.  Save. */
	    int64_t which = stinger_int64_fetch_add (&nslot, 1);
	    if (which < ninsert_remaining) {        /* Can be racy. */
	      has_slot[which] = tmp;
	      kslot[which] = k;
	    }
	  }
	}

	if (nslot < ninsert_remaining) {
	/* Gather any trailing slots. */
	MTA ("mta assert nodep")
	  for (; endk < STINGER_EDGEBLOCKSIZE; ++endk) {
	    int64_t which = stinger_int64_fetch_add (&nslot, 1);
	    if (which < ninsert_remaining) {        /* Can be racy. */
	      has_slot[which] = tmp;
	      kslot[which] = endk;
	    }
	  }
      }
    }
  }

  /* At this point, know how many need inserted & how many slots
     already are available. */

  if (ninsert_remaining > nslot) {
    /* Need more edge blocks. */
    int64_t neb =
      (ninsert_remaining + STINGER_EDGEBLOCKSIZE - 1) / STINGER_EDGEBLOCKSIZE;
    eb_index_t *ebs;
    ebs = xmalloc (neb * sizeof (*ebs));

    new_ebs (ebs, neb, type, from);
    for (int64_t kb = 0; kb < neb - 1; ++kb)
      ebpool[ebs[kb]].next = kb + 1;
    ebpool[ebs[neb - 1]].next = 0;
    *prev_loc = ebs[0];
    push_ebs (G, neb, ebs);

    MTA ("mta assert nodep")
      for (int64_t kb = 0; kb < neb; ++kb) {
        struct stinger_eb *eb = ebpool + ebs[kb];
        MTA ("mta assert nodep") MTASTREAMS ()
          for (int64_t k = 0;
               nslot < ninsert_remaining
                 && k <
                 STINGER_EDGEBLOCKSIZE; ++k) {
            int64_t which = stinger_int64_fetch_add (&nslot, 1);
            if (which < ninsert_remaining) {        /* Can be racy. */
              has_slot[which] = eb;
              kslot[which] = k;
            }
          }
      }

    free (ebs);
  }

  if (ninsert_remaining) {
    nslot = 0;
    for (int64_t k = 0; k < ninsert; ++k) {
      const int64_t w = insert[k];
      if (w >= 0) {
        int64_t ks;
        struct stinger_eb *restrict eb;
        int64_t which = stinger_int64_fetch_add (&nslot, 1);

        assert (which < ninsert_remaining);
        ks = kslot[which];
        eb = has_slot[which];

        assert (ks < STINGER_EDGEBLOCKSIZE);
        assert (ks >= 0);
        assert (eb);
        /* Breaking atomicity => assert may break. */
        update_edge_data (G, eb, ks, w,
                          (weight ? weight[k] : 1), timestamp);
      }
    }
  }

  free (kslot);
  free (has_slot);
  return (nrem + ninsert_remaining) > 0;
}
MTA("mta parallel default")

/** @brief DEPRECATED Process edge removals and insertions as a batch.
 *
 *  DEPRECATED
 *  Takes its input from stinger_sort_actions().  Takes a sorted batch of edge
 *  insertions and removals and processes them in parallel in the data structure.
 *
 *  @param G The STINGER data structure
 *  @param type Edge type
 *  @param timestamp The current timestamp for this batch
 *  @param n Number of incident vertices in the batch
 *  @param insoff For each incident vertex, the offset into the actions array of insertions
 *  @param deloff For each incident vertex, the offset into the actions array of deletions
 *  @param act The sorted actions array
 *  @return The number of incident vertices
 */
int64_t
stinger_remove_and_insert_batch (struct stinger * G, int64_t type,
                                 int64_t timestamp, int64_t n,
                                 int64_t * insoff, int64_t * deloff,
                                 int64_t * act)
{
  /* Separate each vertex's batch into deletions and insertions */
  OMP ("omp parallel for")
  MTA ("mta assert nodep")
  for (int k = 0; k < n; k++) {
    const int64_t i = act[2 * deloff[k]];
    const int64_t nrem = insoff[k] - deloff[k];
    const int64_t nins = deloff[k + 1] - insoff[k];
    int64_t *restrict torem, *restrict toins;

    torem = (nrem + nins ? xmalloc ((nrem + nins) * sizeof (*torem)) : NULL);
    toins = (torem ? &torem[nrem] : NULL);

    MTA ("mta assert nodep")
    for (int k2 = 0; k2 < nrem; ++k2) {
      torem[k2] = act[2 * (deloff[k] + k2) + 1];
      assert (act[2 * (deloff[k] + k2)] == i);
    }

    MTA ("mta assert nodep")
    for (int k2 = 0; k2 < nins; k2++) {
      toins[k2] = act[2 * (insoff[k] + k2) + 1];
      assert (act[2 * (insoff[k] + k2)] == i);
    }

/* Do the update in the data structure */
    stinger_remove_and_insert_edges (G, type, i, nrem, torem, nins, toins,
                                       NULL, timestamp);
      if (torem != NULL)
	free (torem);
    }

  return n;
}

/** @brief DEPRECATED Copy typed adjacencies of a vertex into a buffer
 *
 *  DEPRECATED
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
void
stinger_gather_typed_successors_serial (const struct stinger *G, int64_t type,
                                        int64_t v, size_t * outlen,
                                        int64_t * out, size_t max_outlen)
{
  size_t kout = 0;

  assert (G);

  if (v < 0) {
    *outlen = 0;
    return;
  }

  STINGER_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(G, v, type) {
    if(kout < max_outlen)
      out[kout++] = STINGER_EDGE_DEST;
  } STINGER_FORALL_EDGES_OF_TYPE_OF_VTX_END();

  *outlen = kout;               /* May be longer than max_outlen. */
}


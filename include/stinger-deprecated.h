#ifndef STINGER_DEPRECATED_H
#define STINGER_DEPRECATED_H

/* These functions are included for backwards compatability, but are not recommended
 * for use.  No guarantees are made about their performance, and they are no longer
 * part of the STINGER API.
 */

int
stinger_remove_and_insert_edges (struct stinger *G,
                                 int64_t type, int64_t from,
                                 int64_t nremove, int64_t * remove,
                                 int64_t ninsert, int64_t * insert,
                                 int64_t * weight, int64_t timestamp);

int64_t
stinger_remove_and_insert_batch (struct stinger * G, int64_t type,
                                 int64_t timestamp, int64_t n,
                                 int64_t * insoff, int64_t * deloff,
                                 int64_t * act);
void
stinger_gather_typed_successors_serial (const struct stinger *G, int64_t type,
                                        int64_t v, size_t * outlen,
                                        int64_t * out, size_t max_outlen);

#endif  /*STINGER-DEPRECATED_H*/


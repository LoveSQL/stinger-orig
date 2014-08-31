/**
 * @mainpage STINGER Documentation
 *
 * @section intro-sec Introduction
 * <a href="http://www.cc.gatech.edu/stinger">STINGER</a> documentation is maintained using 
 * <a href="http://www.doxygen.org">Doxygen</a>.
 *
 * @section module-sec Module overview
 * Everything you need to know is contained in a few header files: 
 * - include/stinger.h:
 *     - STINGER Core API
 *     - Create and delete STINGERs
 *     - Insert, remove, and increment edges
 *     - Get and set metadata on vertices, edges
 *     - Get counts of vertices and edges
 *     - Get adjacencies (successors and predecessors) of vertices
 *     - Save / load STINGER to / from disk
 *     - C Macros to aid in fast graph traversal
 * - include/stinger-utils.h
 *     - STINGER utilities
 *     - Read gen-streams formatted data
 *     - Batch insert data into STINGER
 *     - sorted / unsorted CSR to STINGER
 *     - STINGER to sorted / unsorted CSR
 *     - Edge lists (simple arrays of edges) to STINGER
 *     - Sorting and comparison utilities
 * - include/stinger-physmap.h
 *     - STINGER physical mapper
 *     - External structure to map arbitrary strings and binary data into STINGER vertex space
 * - include/stinger-iterator.h
 *     - STINGER filtering iterator
 *     - Serially iterates the set of edges matching the user-specified filter
 *     - Filter on any combination of
 *         - Set of vertices
 *         - Set of vertex types
 *         - Set of edge types
 *         - First and recent timestamps
 *         - User-defined filter
 *
 * Extra utilities can be found here:
 * - include/timer.h
 *     - System Timers
 * - include/xmalloc.h
 *     - Memory Allocation Wrappers
 */

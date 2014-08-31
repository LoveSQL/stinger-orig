<meta http-equiv='Content-Type' content='text/html; charset=utf-8'>

STINGER: Getting Started
========================

Welcome.  This document was written corresponding to the R633 stable release of STINGER.  It is intended to be 
a short but reasonably comprehensive guide to the code base, including how to build the included examples and run
a few basic demos, how to compile STINGER as a static library, and how to get started using STINGER in a simple project.

Target Systems
--------------

This version of STINGER is designed to run on large memory Intel-compatible x86 server platforms and the 
Cray XMT and is capable of processing millions to billions of vertices and edges on these platforms.  It is also capable
of running on your laptop or desktop on graphs of millions of edges and vertices; however, you may need to
adjust the data structure's maximum storage size to fit the memory of your machine. 

To adjust STINGER's allocation size, you can tune the maximum number of vertices in the graph, the number of 
edges in a STINGER block (recommended to be between 12 and 32), and the maximum number of edge types
(recommended to between 0 and 15 or so due to the storage required for edge type indices) in include/stinger-defs.h:

    /** Maximum number of vertices in STINGER */
    #if !defined(STINGER_MAX_LVERTICES)
      #if defined (__MTA__)
        /* XXX: Assume only 2**25 vertices */
        #define STINGER_MAX_LVERTICES (1L<<27) /* 2 ^ 27 - 134 million*/
      #else
        /* much smaller for quick testing */
        #define STINGER_MAX_LVERTICES (1L<<21) /* 2 ^ 21 - 2 million */
      #endif
    #endif

    /** Number of edges in each edge block */
    #define STINGER_EDGEBLOCKSIZE 32
    /** Number of edge types */
    #define STINGER_NUMETYPES 1

    
The total maximum number of edges is determined by the maximum number of edge blocks which is defined
in include/stinger-internal.h as:

    #define EBPOOL_SIZE (STINGER_MAX_LVERTICES*4)

**In the event you see the message "XXX: eb pool exhausted",** you should adjust these numbers to create
more edges (increase STINGER\_MAX\_LVERTICES or increase the factor in EBPOOL\_SIZE).

**In the event you see a message like "src/stinger.c:786: stinger\_insert\_edge: Assertion `from < (1L<<21)' failed.",** 
you should increase STINGER\_MAX\_LVERTICES.

Keep in mind that increasing these numbers could cause STINGER exceed the limits of your computer's RAM.

Directory Structure
-------------------

    .
    ├── compat - XMT compatability libraries
    ├── doc - basic README, LICENSE, Doxygen configuration
    ├── examples - contains various standalone executable examples using STINGER
    │   ├── breadth_first_search
    │   ├── clustering_coefficients_streaming
    │   ├── community-reagglomeration
    │   ├── connected_components_static
    │   ├── connected_components_streaming
    │   ├── EMPTY_PROJECT_START_HERE
    │   ├── insert_remove_benchmark
    │   ├── multi_contract_clustering
    │   └── new_streaming_clustering_coefficients
    ├── genstreams - R-MAT graph and stream generator
    ├── include - all header files for STINGER
    ├── java
    │   ├── doxy
    │   └── stinger
    ├── obj
    ├── python
    │   └── doxy
    ├── server
    ├── src
    └── test


Building STINGER
----------------

The first step in building STINGER is to select a compiler appropriate for your platform.
STINGER is built using GNU-compatible Makefiles with support for gcc and clang on x86 \*nix, OS X, and Cygwin, and 
the Cray compiler for the XMT.
Include files are given for four configurations:

1. make.inc-crayxmt
2. make.inc-gcc
3. make.inc-gcc-openmp
4. make.inc-clang

Create a link called "make.inc" to your desired configuration.

    ln make.inc-gcc-openmp make.inc

Then call "make" to build all of the default targets.

    make

Building a static library
-------------------------

Call make on the lib target to produce libstinger.a.  Compile this into your programs and use the headers
in the include directory.

    make lib

To use the produced libstinger.a file, you will need to link against it when compiling/linking your 
code.  For example, if you have a file main.c that uses STINGER, this will build your application:

    gcc -std=c99 -Iinclude main.c -L. -lstinger -lrt

Note also that if you built libstinger.a using OpenMP, you will need the OpenMP flag:

    gcc -std=c99 -Iinclude main.c -L. -lstinger -lrt -fopenmp

On Linux, the -lrt flag is required, but it may be optional on other platforms.

Building Java and Python Libraries
----------------------------------

Requirements:

* Relevant library and headers (JDK and/or Python)
* [SWIG](www.swig.org)
* [Doxygen (for documentation)](www.doxygen.org)
* Select configuration as in Building STINGER

STINGER has unsupported toy interfaces to Python and Java.
These interfaces are configured and constructed using the SWIG automated tool.  
As such, they are for convenience, and not necessarily for performance.  The paths to the 
headers of your libraries should be configured in 
"Makefile" in the variables JAVA\_LIB and/or PYTHON\_LIB.

Keep in mind that when using these interfaces, you will need to use the SWIG array classes
in your respective language to pass in arrays that would otherwise be passed as pointers in C.

### Java

To build the java library, call make on the java target:

    make java

When the build completes, in the java directory you should have a libstinger.so library that implements the 
Java Native Interface for STINGER.  There is an example of how to use this in java/main.java.  To run this example
from the java directory:

    javac main.java
    java main

### Python

To build the Python library, call make on the python target:

    make python

When the build completes, in the python directory you should have a \_stinger.so library that implements the 
Python native module for STINGER.  There is an example of how to use this in python/main.py.  To run this example
from the python directory:

    python main.py

Generating Data
---------------

Most included examples are designed to operate on synthetic graph data in a compressed sparse row binary format and stream
data in a binary tuple format (see snarf\_graph and snarf\_actions in src/stinger-utils.c).  STINGER comes with
two tools to produce these graphs, gen-streams and rmatter.  Both are built by the default make target. 
The most commonly used is gen-streams.  

### gen-streams
Gen-streams produces a synthetic graph using an implementation of the Recursive MATrix (R-MAT) algorithm and also
produces a stream of edge insertions and deletions using the same algorithm.  Deletions are sampled from previously
generated edges.  The parameters for gen-streams include scale (--scale=20) which specifies the number of vertices 
in the output graph as 2^SCALE, edgefactor (--edgefact=8) which specifies the average degree of the vertices in the initial
graph, the number of actions (--nact=1000000) which specifies how many updates will be in the output update stream, and many 
other options which can be seen through the -? flag.  

An example of running gen-streams to create a graph with one million
vertices and 16 million undirected edges as well as an actions stream containing one million actions:

    ./gen-streams --scale=20 --edgefact=16 --nact=1000000 graphoutfile.bin actionsoutfile.bin

Note that the output filenames and all other parameters are optional.


### rmatter
Rmatter is also an R-MAT graph generator similar to gen-streams.  It is only compiled on x86.  rmatter's random
number generation uses drand48() and is less mathematically sound than the approach in gen-streams; however, 
rmatter is parallel, works well for large graphs, and generally runs faster.  To produce the same size of graph and 
action stream:

    ./rmatter -s 20 -e 16 -n 1000000 -g graphoutfile.bin -a actionsoutfile.bin

Running Examples
----------------

### Default: Static Connected Components
The executable that is built from main.c in the trunk is a simple example of how basic algorithms can
be implemented on graphs STINGER.  In particular, it loads a static graph into STINGER, labels the connected components
using sequential and parallel versions of implementations based on the Shiloach-Vishkin algorithm, and prints out the 
times and sizes.  This example is built by the default target.

    ./main graphoutfile.bin actionsoutfile.bin

### Streaming clustering coefficients
There are several example implementations of streaming clustering coefficients included in this distribution
under examples.  The default streaming clustering coefficients target is in examples/clustering\_coefficients\_streaming.
To run it, simply run make or build the make target streaming-clustering-coefficients, then call 
streaming\_clustering\_coefficients on the input graph and actions stream:

    ./streaming_clustering_coefficients graphoutfile.bin actionsoutfile.bin

This will calculate the initial local clustering coefficients of the graph, then update the clustering coefficients 
while inserting and removing batches of updates.  You may specify the number of updates and the size of each update
using the parameters -n and -b, respectively:

    ./streaming_clustering_coefficients -n 10 -b 100000 graphoutfile.bin actionsoutfile.bin

There is a slightly different version of the code in examples/new\_streaming\_clustering\_coefficients that has been 
more rigorously tested.

### Community Detection / Clustering
There are examples of clustering algorithms that ship in this distribution.  The first is a streaming community detection
implementing agglomerative clustering with re-agglomeration to handle graph updates. It can be found in 
examples/community-reagglomeration.  To build it, simply run make, then run it with similar parameters to
the clustering coefficients example above:

    cd examples/community-reagglomeration
    make
    ./main -n 10 -b 100000 ../../graphoutfile.bin ../../actionsoutfile.bin

The other clustering is a static in-place agglomerative approach using multi-level tree contractions. 
It is located in examples/multi\_contract\_clustering and can be built and run similarly.  It will print
out the mapping of vertices to clusters upon completion and contains some example code of how to load
the text format used in the most recent DIMACS graphs. 

### Streaming Connected Components
The examples/connected\_components\_streaming directory contains an example of how to track and 
relabel the components of a vertex using a triangle counting heuristic to account for deletions.  The code can
be compiled and run similar to the previous examples.

### Insert / Remove Benchmark
This program simple loads the graph and actions stream, builds the STINGER graph, and then begins inserting and 
removing batches of edges from the update stream as quickly as possible.  On high-end x86 platforms and the Cray XMT,
this benchmark has been able to reach over 3 million updates per second on graphs with billions of edges.
It is built and run similar to the previous examples.

    cd examples/insert\_remove\_benchmark
    make
    ./main -n 10 -b 100000 ../../graphoutfile.bin ../../actionsoutfile.bin

Programming For STINGER
-----------------------

The key functions for STINGER are defined in include/stinger.h. This includes creating and deleting a STINGER, inserting, 
incrementing, and removing edges, getting and setting the properties of edges and vertices, functions for copying data into
and out of STINGER, and more.  The utilities file (include/stinger-utils.h) contains import and export functions for CSR and 
edge lists, reading functions for STINGER binary data, and other utility functions. Wrappers for atomic functions on your 
platform can be found in include/stinger-atomics.h, and lastly, x86 emulations of the full-empty bit semantics available
on the XMT can be found in include/x86-full-empty.h.

### Creating and Deleting

Creating an empty STINGER is as easy as calling stinger\_new() from stinger.h:

    struct stinger * S = stinger_new();

This creates an empty STINGER with a fixed maximum size in terms of both edges and vertices (see
previous section Target Systems for information on how and where this is set).

To destroy the STINGER and free its memory, call stingeri\_free\_all() from stinger.h:

    S = stinger_free_all(S);

### Inserting, Incrementing, and Removing Edges

STINGER stores all edges in a directed manner; however, to the author's knowledge, all work
done on STINGER so far assumes that the graph is indirected and that each directed edge will 
have a corresponding edge in the other direction.  Functions for manipulating directed edges are 
mentioned below.

Edges in STINGER have a source vertex, a destination vertex, edge type, edge weight, and 
two timestamps.  The type is an integer between 0 and STINGER\_NUMETYPES as defined in stinger-defs.h.
The weight and timestamps are 64-bit integers.  The weight and its meaning are up to the user.  The 
timestamps are also user determined according to certain rules; however, they cannot be INT64\_MAX or
INT64\_MIN as these are reserved for internal use.  The default behavior for the two timestamps
is similar to the creation and modification times in a file system.  When inserting or incrementing an edge, 
if the edge did not previously exist, both the time\_first and time\_recent timestamps will be set to the given timestamp.
If the edge previously existed, the recent timestamp will be updated if the given time is greater.

Edge insertions will overwright any existing edge weights.  Increments will create an edge and set its weight if it didn't 
previously exist.  Otherwise, increments will atomically add the given weight to the existing weight.

Insertions use the function stinger\_insert\_edge().

    /* insert random edges of type 0 */
    for(int i = 0; i < 1000; i++) {
      int64_t source = rand();
      int64_t dest = rand();
      int64_t weight = 1;

      stinger_insert_edge(S, 0, source, dest, weight, time(NULL));
    }

Increments use the same format, but use the function stinger\_incr\_edge()

    stinger_incr_edge(S, 0, source, dest, weight, time(NULL));

Removing an edge will delete the edge completely from the data structure.  This function does not accept a weight
or time.

    stinger_remove_edge(S, 0, source, dest);

If you would like to update the edge timestamp without chaning the weight, use this function:

    stinger_edge_touch (S, source, dest, type, time);

#### Undirected Edge Functions
**All three of these functions manipulate directed edges,** but their are corresponding edge pair
functions that will perform the insertion, increment, or removal operation in both directions for
undirected relationships.  These follow the same conventions as the directed functions above:

    stinger_insert_edge_pair(S, 0, source, dest, weight, time(NULL));
    stinger_incr_edge_pair(S, 0, source, dest, weight, time(NULL));
    stinger_remove_edge_pair(S, 0, source, dest);

### Getting Properties of Edges

There are convenience functions for getting the properties of individual edges, but it is recommended 
that you aggregate examining the properties through traversal or some other method if performance is 
as concern unless absolutely necessary as each lookup can have a cost of _O(deg(_source_))_.

    stinger_edgeweight(S, source, dest, type);

    stinger_edge_timestamp_first(S, source, dest, type);

    stinger_edge_timestamp_recent(S, source, dest, type);

### Adding and Removing Edges in Bulk

All of the insert, increment, and remove functions and the pair variants are thread safe.  If there
are already edges in your graph, the fastest way to insert and remove edges is to call these functions
in a parallel loop or other parallel context.  The operation is complete when the function returns 
(keeping in mind that STINGER's consistency model is fairly loose).  If your graph is empty,
and your data is in an easily digestible format (such as CSR or an edge list), there are functions
to fill an empty STINGER from these inputs.

For example, from compressed sparse row (CSR), simply call stinger\_set\_initial\_edges() as below, 
where NV is the number of
vertices in your CSR (off should be nv + 1 in length), to is the destination vertex packed array, weight is
the weight of the corresponding edge to the destination vertex, time and first\_time are similar.  If
weight is NULL, weights will be set to 1.  If one time is NULL, the other array is used for both timestamps,
but if both time arrays are NULL, the single\_ts is used for all timestamps.  Note that this function
creates all edges with the same given type.

    stinger_set_initial_edges (S, nv, etype, off, to, weight, time, first_time, single_ts);

If your edge is in an edge list form (source vertices, destination vertices, weights, timestamps, etc. all in
separate arrays), you can skip calling STINGER new and build a STINGER directly from the edge\_list\_to\_stinger()
function call in stinger-utils.h:

    struct stinger * S = edge_list_to_stinger (nv, ne, source, dest, weight, timerecent, timefirst, single_ts); 

The weight, timerecent, and timefirst follow the same conventions.  All edges will be assumed to be of
type zero.


### Setting and Getting Vertex Metadata

At the global scale, STINGER can tell you how many vertices are active (must have at least one
outgoing edge) and which active vertex has the highest vertex ID:

    stinger_max_active_vertex(S);

    stinger_num_active_vertices(S);

This can be useful when your algorithm needs to iterate over the vertices from 0 to NV (instead of
using STINGER\_MAX\_LVERTICES).

There are functions for accessing and setting the type and weight metadata associated with each vertex.
There is no maximum on the number of vertex types (as these are unindexed), and the meaning and use of
both the types and weights of vertices are completely up to the user.

    stinger_vweight (S, vertexID);
    stinger_vtype (S, vertexID);
    stinger_set_vweight (S, vertexID, weight);
    stinger_set_vtype (S, vertexID, type);

STINGER tracks the in degree and out degree of each vertex, so the provided functions simply perform
an _O(_1_)_ lookup for these values.  The typed out degree function on the other hand has a cost of
_O(deg(_v_))_.

    stinger_outdegree (S, vertexID);
    stinger_indegree (S, vertexID);
    stinger_typed_outdegree (S, vertexID, type);

### Copying out Data

STINGER provides functions to copy the adjacency list of a vertex out into separate packed arrays via
stinger\_gather\_successors().  The weight, timestamps, and types outputs are optional and can be NULL.
These can be filtered to only a specific type via the stinger\_gather\_typed\_successors() function, but 
this function only provides the vertexIDs.  Both functions have a cost of _O(deg(_v_))_. 
There is also a gather typed predecessors function, but this has a cost of _O(_E_)_.

    stinger_gather_successors(S, source, outlength, dest, weight, timefirst, timerecent, type, max_outlen);
    stinger_typed_gather_successors(S, type, source, outlength, dest, max_outlen);
    stinger_gather_typed_predecessors(S, type, source, outlength, dest, max_outlen);

Similar to these functions, there is a has typed successor function which checks for an out edge of 
a given type to a specific vertex with cost _O(deg(_v_))_.

    stinger_has_typed_successor (S, type, source, dest);

STINGER can export its internal data to CSR format in both sorted and unsorted form.  Weight time and type
fields are all optional and may be NULL.  Note that STINGER will allocate the output arrays, but it is
the user's responsibility to free each of them.

    stinger_to_sorted_csr (S, nv, out, ind, weight, timefirst, timerecent, type);
    stinger_to_unsorted_csr (S, nv, out, ind, weight, timefirst, timerecent, type);

### Saving and Restoring a STINGER

STINGER provides checkpoint functions to save the STINGER structure and all data contained in it to disk
and to reload that data.

    stinger_save_to_file (S, max_active_vtx+1, "/path/to/file");
    stinger_open_from_file ("/path/to/file", &nv_out, &S_out);

### Traversal

STINGER provides convenient traversal macros for abstracting away the details of accessing the data 
structure to make writing algorithms clearer.  These algorithms focus on iterating over edges, either
all edges in the graph of a specific type, or all edges of a vertx (optionally of a specific type).

    STINGER_FORALL_EDGES_OF_VTX_BEGIN(STINGER_,VTX_) {
    } STINGER_FORALL_EDGES_OF_VTX_END();

    STINGER_PARALLEL_FORALL_EDGES_OF_VTX_BEGIN(STINGER_,VTX_) {
    } STINGER_PARALLEL_FORALL_EDGES_OF_VTX_END();

    STINGER_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(STINGER_,TYPE_,VTX_) {
    } STINGER_FORALL_EDGES_OF_TYPE_OF_VTX_END();

    STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_BEGIN(STINGER_,TYPE_,VTX_) {
    } STINGER_PARALLEL_FORALL_EDGES_OF_TYPE_OF_VTX_END();

    STINGER_FORALL_EDGES_BEGIN(STINGER_,TYPE_) {
    } STINGER_FORALL_EDGES_END();

    STINGER_PARALLEL_FORALL_EDGES_BEGIN(STINGER_,TYPE_) {
    } STINGER_PARALLEL_FORALL_EDGES_END();

As you can see above, each traversal macro also provides its own parallel variant.  For parallel code,
we recommend using the platform-independent atomic operations in stinger-atomics.h and the full-empty
semantics in x86-full-empty.h.

When using one of these macro pairs, the code inside should be written from the perspective of one
edge at a time.  There are macros to access each property of an edge from inside the loop:

    STINGER_EDGE_SOURCE
    STINGER_EDGE_DEST
    STINGER_EDGE_TYPE
    STINGER_EDGE_WEIGHT
    STINGER_EDGE_TIME_FIRST
    STINGER_EDGE_TIME_RECENT

For example, to add up all of the weights of edges of vertex 5 that were created before time 7 and
updated after time 9:

    int64_t total = 0;

    STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, 5) {
      if(STINGER_EDGE_TIME_FIRST < 7 && STINGER_EDGE_TIME_RECENT > 9) {
        total += STINGER_EDGE_WEIGHT;
      }
    } STINGER_FORALL_EDGES_OF_VTX_END();

Similarly, to find the lowest total scoring edge of type 0 edges in a graph:

    /* assuming we are given some 0 to NV-1 array of double scores */

    double min_score = DBL_MAX;
    STINGER_PARALLEL_FORALL_EDGES_BEGIN(S, 0) {
      double score = scores[STINGER_EDGE_SOURCE] + scores[STINGER_EDGE_DEST];
      if(score < min_score) {
        OMP("omp critical")
        {
            if(score < min_score)
              min_score = score;
        }
      }
    } STINGER_PARALLEL_FORALL_EDGES_END();

Note here the protection applied to where the score is changed for the parallel case.  This critical
section will be compiled out if OpenMP is not enabled.

### Skeleton Code
The default distribution of STINGER contains a skeleton code in examples/EMPTY\_PROJECT\_START\_HERE.
This code handles parsing command line parameters, loading an initial static graph from the STINGER binary format, and inserting
and removing the edge updates.  Calling make and ./main with the same parameters as above from inside of this directory
will build and run the skeleton code.

How To Build Documentation
--------------------------

Most STINGER functions intended for external use are documented with Doxygen-style comments.  You can generate the 
documentation by installing [doxygen](http://www.doxygen.org), then calling doxygen in the main directory.
This creates html pages in the doc/html directory.  See index.html.


How To Get Help
---------------

This implementation of STINGER has been developed and supported by a handful of individuals in David Bader's
High Performance Computing group at Georgia Tech.  You can get support through http://www.cc.gatech.edu/~bader,
www.stingergraph.com, or through emailing the developers. 

David Ediger dediger@gatech.edu
Rob McColl rmccoll3@gatech.edu
Jason Riedy jason.riedy@cc.gatech.edu

Licensing
---------

STINGER is licensed under the 3-clause BSD license (see doc/LICENSE), which means that you may redistribute it
as you see fit, modify it to your heart's content, and in general do whatever you want to with the code 
base without having to tell anyone about it, without being required to release your subsequent code, and without
any obligation.  

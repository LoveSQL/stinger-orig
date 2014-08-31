#ifndef GRAPH_WITHOUT_TYPES_H
#define GRAPH_WITHOUT_TYPES_H
#include "stinger-tests.h"
#include <ctime>

#define MIN_NV 10
#define MAX_NV 20

// TODO: empty graph test, remove, remove and insert, reinsert edge, timestamp test
// TODO: iterate through edges no matter what there type is




#define GRAPH_SETUP_CLEANUP_METHODS(graph_var) \
void graph_var##_setup() {\
	graph = stinger_new();\
	NV = graph_var##_NV;\
	NE = graph_var##_NE;\
	graph_edges = &graph_var##_edges;\
	graph_source = graph_var##_source;\
	graph_dest = graph_var##_dest;\
	graph_types = graph_var##_types;\
	graph_wgts = graph_var##_wgts;\
	graph_timestamps = graph_var##_timestamps;\
\
	for (int i = 0; i < graph_var##_NE; ++i) {\
\
		graph_var##_edges.insert(\
			pair<edge_t, edge_info_t> (\
			edge_t(graph_var##_source[i],graph_var##_dest[i]),\
			edge_info_t(graph_var##_types[i], graph_var##_wgts[i], graph_var##_timestamps[i]) )\
			);\
	}\
}\
\
void graph_var##_cleanup()\
{\
	stinger_free(graph);\
	graph_var##_edges.clear();\
}\

const int64_t graph0_NV = 30;
const int64_t graph0_NE = 49;
int64_t graph0_source[] = {
	0,  0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,
	3,  3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,
	5,  5,  5,  5,  5,  6,  7,  8,  9,  9,  10, 11,
	12, 13, 17, 18, 18, 18, 19, 20, 22, 23, 24, 26,
	27
};
int64_t graph0_dest[] = {
	1,  3,  2,  4,  7,  8,  10, 11, 15, 3,  5,  8,
	4,  6,  9,  12, 14, 15, 16, 18, 21, 5,  8,  10,
	12, 15, 17, 19, 24, 8,  12, 10, 16, 24, 13, 15,
    15, 19, 22, 23, 26, 27, 28, 29, 23, 28, 29, 29,
	28
};
int64_t graph0_types[] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0};
int64_t graph0_wgts[] = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
	 		 	  1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1,
				  1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1
				 };
int64_t graph0_timestamps[] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0};

edgemap_t graph0_edges;

const int64_t graph1_NV = 30;
const int64_t graph1_NE = 49;
int64_t graph1_source[] = {
	0,  0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,
	3,  3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,
	5,  5,  5,  5,  5,  6,  7,  8,  9,  9,  10, 11,
    12, 13, 17, 18, 18, 18, 19, 20, 22, 23, 24, 26,
	27
};
int64_t graph1_dest[] = {
	1,  3,  2,  4,  7,  8,  10, 11, 15, 3,  5,  8,
	4,  6,  9,  12, 14, 15, 16, 18, 21, 5,  8,  10,
	12, 15, 17, 19, 24, 8,  12, 10, 16, 24, 13, 15,
    15, 19, 22, 23, 26, 27, 28, 29, 23, 28, 29, 29,
	28
};
int64_t graph1_types[] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0};
int64_t graph1_wgts[] = {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16,
	 		 	  17,18,19,20, 101,102,103,104, 105,106,107,108, 109,1010,1011,1012,
				  1013,1014,1015,1016, 1017,1018,1019,1020, 1021,1022,1023,1024, 1025,1026,1027,1028, 1029
				 };
int64_t graph1_timestamps[] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0};

edgemap_t graph1_edges;

const int64_t empty_graph_NV = 0;
const int64_t empty_graph_NE = 0;
int64_t empty_graph_source[] = {};
int64_t empty_graph_dest[] = {};
int64_t empty_graph_types[] = {};
int64_t empty_graph_wgts[] = {};
int64_t empty_graph_timestamps[] = {};

edgemap_t empty_graph_edges;


GRAPH_SETUP_CLEANUP_METHODS(graph0)
GRAPH_SETUP_CLEANUP_METHODS(graph1);
GRAPH_SETUP_CLEANUP_METHODS(empty_graph);

int64_t graph_random_NV;
int64_t graph_random_NE;
int64_t* graph_random_source;
int64_t* graph_random_dest;
int64_t* graph_random_types;
int64_t* graph_random_wgts;
int64_t* graph_random_timestamps;

edgemap_t graph_random_edges;

void graph_random_setup() {
        //srand(time(0));
	graph = stinger_new();

	graph_random_NV = MIN_NV + rand() % (MAX_NV-MIN_NV);
	graph_random_NE = graph_random_NV + rand() % (graph_random_NV*graph_random_NV - graph_random_NV);

	graph_random_source = new int64_t[graph_random_NE];
	graph_random_dest = new int64_t[graph_random_NE];
	graph_random_types = new int64_t[graph_random_NE];
	graph_random_wgts = new int64_t[graph_random_NE];
	graph_random_timestamps = new int64_t[graph_random_NE];

	graph_edges = &graph_random_edges;
        int64_t nignores = 0;
	for (int64_t i = 0; i < graph_random_NE; ++i) {
		int iter = 0;
		const int maxiter = 100;
		do {
			graph_random_source[i-nignores] = rand() % (graph_random_NV);
			graph_random_dest[i-nignores] = rand() % (graph_random_NV);
			if(graph_random_source[i-nignores] > graph_random_dest[i-nignores]){
			  int64_t tmp = graph_random_source[i-nignores];
			  graph_random_source[i-nignores] = graph_random_dest[i-nignores];
			  graph_random_dest[i-nignores] = tmp;
			}
			++iter;
		} while (graph_random_edges.find(edge_t(graph_random_source[i-nignores],graph_random_dest[i-nignores])) != graph_random_edges.end() && iter < maxiter);
		if (iter >= maxiter) {
		  nignores++;
	          //srand(time(0));
		  continue;
		}

		graph_random_types[i-nignores] = 0;
		graph_random_wgts[i-nignores] = rand() % graph_random_NV;
		graph_random_timestamps[i-nignores] = 1; //rand() % graph_random_NV;

		graph_random_edges.insert(
			pair<edge_t, edge_info_t> (
			edge_t(graph_random_source[i-nignores],graph_random_dest[i-nignores]),
			edge_info_t(graph_random_types[i-nignores], graph_random_wgts[i-nignores], graph_random_timestamps[i-nignores]) )
			);
	}


	NV = graph_random_NV;
	NE = graph_random_NE - nignores;
	graph_source = graph_random_source;
	graph_dest = graph_random_dest;
	graph_types = graph_random_types;
	graph_wgts = graph_random_wgts;
	graph_timestamps = graph_random_timestamps;
}

void graph_random_cleanup() {

	delete []graph_random_source;
	delete []graph_random_dest;
	delete []graph_random_types;
	delete []graph_random_wgts;
	delete []graph_random_timestamps;

	stinger_free(graph);
	graph_random_edges.clear();
}

int graph_num_tests_untyped = 3;
void (*graph_setup[])(void) = {graph0_setup, graph1_setup, empty_graph_setup};
void (*graph_cleanup[])(void) = {graph0_cleanup, graph1_cleanup, empty_graph_cleanup};
char* graph_name[] = {"graph0", "graph1", "empty_graph"};

#endif

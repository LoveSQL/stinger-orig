extern "C" {
#include <check.h>
}

#include <cstdlib>
#include <cstring>
#include <set>
#include <map>
#include "stinger-tests.h"
#include "test-register.h"
#include <string>

#include <map>
#include <algorithm>
#include <string>
#include <cassert>
#include <cstdio>

#define edge_t pair<int64_t,int64_t>

#define NUM_ELEMENTS(ARR) (sizeof ((ARR)) / sizeof (int64_t))

#define REGISTER(_TEST_NAME_,_FUNC_NAME_) do {\
  test_names.insert( pair<string, TEST_ID> (string(#_TEST_NAME_) , TEST_##_TEST_NAME_ ));\
  test_ids.insert( pair<TEST_ID, string>(TEST_##_TEST_NAME_,string(#_TEST_NAME_)) );\
  test_funcs[TEST_##_TEST_NAME_] =  test_##_FUNC_NAME_;\
} while (0);


#define CHECK_ALL 0xFFFFFFFF
#define CHECK_TIMESTAMPS 0x000000FF
#define CHECK_WEIGHTS 0x0000FF00


using namespace std;

static map<string, TEST_ID> test_names;
static map<TEST_ID, string> test_ids;
static void (*test_funcs[1000])(int);

typedef map<edge_t, edge_info_t> edgemap_t;

struct stinger* graph;
edgemap_t* graph_edges;
int64_t NV;
int64_t NE;
int64_t *graph_source;
int64_t *graph_dest;
int64_t *graph_types;
int64_t *graph_wgts;
int64_t *graph_timestamps;

char tests_to_run[1000] = { 0 };

#include "graph_without_types.h"

int num_random_tests = 1;

static void init_test_registry();
static TEST_ID find_test_id(char* test_name);

void test_NULL(int a);

int check_edges_exist(int check_flag) {

  edgemap_t graph_edges_copy;
  edgemap_t::iterator iter;
  for (iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
    graph_edges_copy.insert( pair<edge_t,edge_info_t>((*iter).first, (*iter).second) );
  }

  int consistency_result = stinger_consistency_check(graph, NV);
  fail_unless(consistency_result == 0, "Consistency check failed: %d\n", consistency_result);

  STINGER_FORALL_EDGES_BEGIN(graph, 0)

    edgemap_t::iterator iter = 
    graph_edges_copy.find(edge_t(STINGER_EDGE_SOURCE, STINGER_EDGE_DEST));


  fail_unless(iter != graph_edges_copy.end(), "Non-Existent edge (%ld,%ld) returned by stinger", STINGER_EDGE_SOURCE, STINGER_EDGE_DEST);

  graph_edges_copy.erase(iter);
  edge_info_t info = (*iter).second;
  fail_unless(info.type == STINGER_EDGE_TYPE, "Edge (%ld,%ld) of invalid type %ld returned, expected: %ld",
      STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_TYPE, info.type);
  if (check_flag & CHECK_WEIGHTS)
    fail_unless(info.wgt == STINGER_EDGE_WEIGHT, "Edge (%ld,%ld) Weight %ld is not correct, expected: %ld",
	STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_WEIGHT, info.wgt);
  if (check_flag & CHECK_TIMESTAMPS)
    fail_unless(info.timestamp == STINGER_EDGE_TIME_RECENT, "Edge (%ld,%ld) Timestamp %ld is not correct, expected: %ld",
	STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_TIME_RECENT, info.timestamp);

  STINGER_FORALL_EDGES_END();

  for (edgemap_t::iterator iter = graph_edges_copy.begin(); iter != graph_edges_copy.end(); ++iter) {
    edge_t edge = (*iter).first;
    edge_info_t info = (*iter).second;
    fail_unless(info.type == 0, "Edge (%ld,%ld) of type %ld missing in STINGER data structure", edge.first,
	edge.second, info.type);
  }

  return 0;
}

int check_edge_pairs_exist() {

  edgemap_t graph_edges_copy;
  STINGER_FORALL_EDGES_BEGIN(graph, 0)

    edgemap_t::iterator iter = 
    graph_edges->find(edge_t(STINGER_EDGE_SOURCE, STINGER_EDGE_DEST));

  if (iter == graph_edges->end()) {

    iter = graph_edges->find(edge_t(STINGER_EDGE_DEST, STINGER_EDGE_SOURCE));
    if (iter == graph_edges->end()) {
      iter = graph_edges_copy.find(edge_t(STINGER_EDGE_SOURCE, STINGER_EDGE_DEST));
      fail_unless(iter != graph_edges_copy.end(), "Non-Existent edge (%ld,%ld) returned by stinger", STINGER_EDGE_SOURCE, STINGER_EDGE_DEST);
      graph_edges_copy.erase(iter);
    }
  }
  else {
    graph_edges_copy.insert( pair<edge_t,edge_info_t>(edge_t(STINGER_EDGE_DEST, STINGER_EDGE_SOURCE),(*iter).second) );
    graph_edges->erase(iter);
  }
  edge_info_t info = (*iter).second;
  fail_unless(info.type == STINGER_EDGE_TYPE, "Edge (%ld,%ld) of invalid type %ld returned, expected: %ld",
      STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_TYPE, info.type);
  fail_unless(info.wgt == STINGER_EDGE_WEIGHT, "Edge (%ld,%ld) Weight is not correct %ld returned, expected: %ld",
      STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_WEIGHT, info.wgt);
  fail_unless(info.timestamp == STINGER_EDGE_TIME_RECENT, "Edge (%ld,%ld) Timestamp is not correct %ld returned, expected: %ld",
      STINGER_EDGE_SOURCE, STINGER_EDGE_DEST, STINGER_EDGE_TIME_RECENT,info.timestamp);

  STINGER_FORALL_EDGES_END();

  for (edgemap_t::iterator iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
    edge_t edge = (*iter).first;
    edge_info_t info = (*iter).second;
    fail_unless(info.type == 0, "Edge (%ld,%ld) of type %ld missing in STINGER data structure", edge.first,
	edge.second, info.type);
  }

  return 0;
}

START_TEST(test_stinger_insert_edge)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  check_edges_exist(CHECK_ALL);
}
END_TEST

START_TEST(test_stinger_insert_edge_pair)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge_pair(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }
  check_edge_pairs_exist();
}
END_TEST

START_TEST(test_stinger_remove_edge)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  int consistency_result = stinger_consistency_check(graph, NV);
  fail_unless(consistency_result == 0, "Consistency check failed: %d\n", consistency_result);

  // A bit of Fuzzy Testing. Randomly select edges to remove.
  int num_edges_to_remove = NE / 4;
  set<int> removed_edges;
  for (int i = 0; i < num_edges_to_remove; ++i) {
    int edge_index = rand() % NE;
    while (removed_edges.find(edge_index) != removed_edges.end()) {
      edge_index = rand() % NE;
    }
    removed_edges.insert(edge_index);

    if(graph_source[edge_index] != graph_dest[edge_index]) {
      fail_unless(stinger_remove_edge(graph, 0, graph_source[edge_index], graph_dest[edge_index]) == 1,
	  "Failed to remove edge %ld,%ld", graph_source[edge_index], graph_dest[edge_index]);

      consistency_result = stinger_consistency_check(graph, NV);
      fail_unless(consistency_result == 0, "Consistency check failed: %d %d\n", consistency_result, i);

      fail_unless(graph_edges->erase(edge_t(graph_source[edge_index], graph_dest[edge_index])) != 0, "Internal Error");
    } else {
      i--;
    }
  }

  check_edges_exist(CHECK_ALL);
}
END_TEST

START_TEST(test_stinger_incr_edge_pair)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge_pair(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int i = 0; i < NE; ++i) {

    if (graph_source[i] >= graph_dest[i])
      continue;

    int64_t wgt_incr = rand() % 100;
    int64_t ts_incr = rand() % 100;

    stinger_incr_edge_pair(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	wgt_incr, graph_timestamps[i]+ts_incr);

    edgemap_t::iterator iter = graph_edges->find( edge_t(graph_source[i],graph_dest[i]) );
    fail_unless(iter != graph_edges->end(), "Internal Error");

    (*iter).second.wgt = graph_wgts[i] + wgt_incr;
    (*iter).second.timestamp = graph_timestamps[i] + ts_incr;
  }

  for (int i = 0; i < NE; ++i) {

    if (graph_source[i] < graph_dest[i])
      continue;

    edgemap_t::iterator iter = graph_edges->find( edge_t(graph_source[i],graph_dest[i]) );
    edgemap_t::iterator iter1 = graph_edges->find( edge_t(graph_dest[i],graph_source[i]) );
    fail_unless(iter != graph_edges->end(), "Edge (%ld,%ld) does not exist", graph_source[i], graph_dest[i]);
    fail_unless(iter1 != graph_edges->end(), "Reverse Edge (%ld,%ld) does not exist", graph_dest[i], graph_source[i]);
    (*iter).second.wgt = (*iter1).second.wgt;
    (*iter).second.timestamp = (*iter1).second.timestamp;
  }

  check_edge_pairs_exist();
}
END_TEST

START_TEST(test_stinger_remove_edge_pair)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge_pair(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  // A bit of Fuzzy Testing. Randomly select edges to remove.
  int num_edges_to_remove = NE / 4;
  set<int> removed_edges;
  for (int i = 0; i < num_edges_to_remove; ++i) {
    int edge_index = rand() % NE;
    while (removed_edges.find(edge_index) != removed_edges.end()) {
      edge_index = rand() % NE;
    }
    removed_edges.insert(edge_index);
    if(graph_source[edge_index] != graph_dest[edge_index]) {
      int  rslt = stinger_remove_edge_pair(graph, 0, graph_source[edge_index], graph_dest[edge_index]);
      fail_unless(rslt == 2,
	  "Failed to remove edge %ld,%ld (code %d)", graph_source[edge_index], graph_dest[edge_index], rslt);
      fail_unless(graph_edges->erase(edge_t(graph_source[edge_index], graph_dest[edge_index])) != 0, "Internal Error");
    }
  }

  check_edge_pairs_exist();
}
END_TEST

START_TEST(test_stinger_incr_edge)
{
  for (int i = 0; i < NE; ++i) {

    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int i = 0; i < NE; ++i) {

    int64_t wgt_incr = rand() % 100;
    int64_t ts_incr = rand() % 100;

    stinger_incr_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	wgt_incr, graph_timestamps[i]+ts_incr);

    edgemap_t::iterator iter = graph_edges->find( edge_t(graph_source[i],graph_dest[i]) );
    fail_unless(iter != graph_edges->end(), "Internal Error");

    (*iter).second.wgt = graph_wgts[i] + wgt_incr;
    (*iter).second.timestamp = graph_timestamps[i] + ts_incr;
  }

  check_edges_exist(CHECK_ALL);
}
END_TEST

START_TEST(test_stinger_total_edges)

  int result = 0;
  for (int i = 0; i < NE; ++i) {

    if(graph_source[i] != graph_dest[i]) {
      result += stinger_insert_edge(graph, graph_types[i], 
	  graph_source[i], graph_dest[i],
	  graph_wgts[i], graph_timestamps[i]);
    }
  }

int orig_nedges = stinger_total_edges(graph);
fail_unless(orig_nedges == result, "After Insertion number of edges returned %ld not same as expected %ld", orig_nedges, result);

// A bit of Fuzzy Testing. Randomly select edges to remove.
int num_edges_to_remove = NE / 4;
set<int> removed_edges;
int remove_result = 0;
for (int i = 0; i < num_edges_to_remove; ++i) {
  int edge_index = rand() % NE;
  while (removed_edges.find(edge_index) != removed_edges.end()) {
    edge_index = rand() % NE;
  }
  removed_edges.insert(edge_index);
  if(graph_source[edge_index] != graph_dest[edge_index]) {
    remove_result++;
    fail_unless(stinger_remove_edge(graph, 0, graph_source[edge_index], graph_dest[edge_index]) == 1,
	"Failed to remove edge %ld,%ld", graph_source[edge_index], graph_dest[edge_index]);
    fail_unless(graph_edges->erase(edge_t(graph_source[edge_index], graph_dest[edge_index])) != 0, "Internal Error");
  }
}

int new_nedges = stinger_total_edges(graph);
fail_unless(new_nedges == orig_nedges-remove_result, "After Removal number of edges returned %ld not same as expected %ld", new_nedges, orig_nedges-remove_result);

// Now add back the edges
set<int>::iterator iter;
for (iter = removed_edges.begin(); iter != removed_edges.end(); ++iter) {
  int edge_index = (*iter);
  stinger_insert_edge(graph, graph_types[edge_index], 
      graph_source[edge_index], graph_dest[edge_index],
      graph_wgts[edge_index], graph_timestamps[edge_index]);
}

new_nedges = stinger_total_edges(graph);
fail_unless(new_nedges == orig_nedges, "After Insertion number of edges returned %ld not same as expected %ld", new_nedges, orig_nedges);

END_TEST

START_TEST(test_stinger_remove_and_insert_edges)
{

  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int64_t i = 0; i < NV; ++i) {

    int64_t from = i;
    int64_t num_edges_to_insert = 0;
    int64_t num_edges_to_remove = 0;

    set<int64_t> missing_vertices;
    for (int j = 0; j < NV; ++j) {
      missing_vertices.insert(j);
    }

    int64_t edges_to_remove[MAX_VTX];
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(graph, from)
      if (rand() % 3 == 0) {

	/*				if (STINGER_EDGE_DEST == 80 && STINGER_EDGE_WEIGHT == 41) {
					printf("here\n");
					int a;
					scanf("%d\n", &a);
					}
	 */
	edges_to_remove[num_edges_to_remove++] = STINGER_EDGE_DEST;
	graph_edges->erase( edge_t(STINGER_EDGE_SOURCE, STINGER_EDGE_DEST) );
      }
      else {
	missing_vertices.erase(STINGER_EDGE_DEST);
      }
    STINGER_FORALL_EDGES_OF_VTX_END();

    int64_t edges_to_insert[MAX_VTX];
    for (set<int64_t>::iterator iter = missing_vertices.begin(); iter != missing_vertices.end(); ++iter) {
      if (rand() % 4 == 0) {
	edges_to_insert[num_edges_to_insert++] = *iter;
	graph_edges->insert( pair<edge_t, edge_info_t>(edge_t(from, *iter), edge_info_t(0,1,100)) );
      }
    }

    stinger_remove_and_insert_edges(graph, 0, from, num_edges_to_remove, 
	edges_to_remove, num_edges_to_insert, edges_to_insert, NULL, 100);

    check_edges_exist(CHECK_ALL);
  }

}
END_TEST

START_TEST(test_stinger_remove_and_insert_edges_batch)
{

  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int64_t i = 0; i < NV; ++i) {

    int64_t from = i;
    int64_t num_edges_to_insert = 0;
    int64_t num_edges_to_remove = 0;

    set<int64_t> missing_vertices;
    for (int j = 0; j < NV; ++j) {
      missing_vertices.insert(j);
    }

    int64_t edges_to_remove[MAX_VTX];
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(graph, from)
      if (rand() % 3 == 0) {

	/*				if (STINGER_EDGE_DEST == 80 && STINGER_EDGE_WEIGHT == 41) {
					printf("here\n");
					int a;
					scanf("%d\n", &a);
					}
	 */
	if(STINGER_EDGE_DEST != STINGER_EDGE_SOURCE) {
	  edges_to_remove[num_edges_to_remove++] = STINGER_EDGE_DEST;
	  graph_edges->erase( edge_t(STINGER_EDGE_SOURCE, STINGER_EDGE_DEST) );
	}
      }
      else {
	missing_vertices.erase(STINGER_EDGE_DEST);
      }
    STINGER_FORALL_EDGES_OF_VTX_END();

    int64_t edges_to_insert[MAX_VTX];
    for (set<int64_t>::iterator iter = missing_vertices.begin(); iter != missing_vertices.end(); ++iter) {
      if (rand() % 4 == 0) {
	if(from != *iter) {
	  edges_to_insert[num_edges_to_insert++] = *iter;
	  graph_edges->insert( pair<edge_t, edge_info_t>(edge_t(from, *iter), edge_info_t(0,1,100)) );
	}
      }
    }

    stinger_remove_and_insert_edges(graph, 0, from, num_edges_to_remove, 
	edges_to_remove, num_edges_to_insert, edges_to_insert, NULL, 100);

    check_edges_exist(CHECK_ALL);
  }

}
END_TEST

START_TEST(test_stinger_gather_typed_successors)
{
  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int i = 0; i < NV; ++i) {

    int from = i;
    size_t outlen;
    int64_t out[MAX_VTX];

    stinger_gather_typed_successors(graph, 0, from, &outlen, out, MAX_VTX);
    fail_unless(outlen >= 0, "Length of output array %ld < 0", outlen);

    edgemap_t graph_edges_copy;
    edgemap_t::iterator iter;
    for (iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      graph_edges_copy.insert( pair<edge_t,edge_info_t>((*iter).first, (*iter).second) );
    }

    for (edgemap_t::iterator iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      edge_t edge = (*iter).first;
      if (edge.first != from  || (edge.first == edge.second) || (*iter).second.type != 0) {
	graph_edges_copy.erase(edge);
      }
    }

    for (int j = 0; j < outlen; ++j) {

      edgemap_t::iterator iter = graph_edges_copy.find(edge_t(from, out[j]));
      fail_unless(iter != graph_edges_copy.end(), "Non Existent edge returned");

      edge_t edge = (*iter).first;

      fail_unless(edge.first == from, "Internal Error");
      fail_unless(edge.second == out[j], "Internal Error");

      graph_edges_copy.erase((*iter).first);
    }

    fail_unless(graph_edges_copy.size() == 0, "Edges Remaining");
  }

}
END_TEST

START_TEST(test_stinger_gather_successors)
{
  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int i = 0; i < NV; ++i) {

    int from = i;
    size_t outlen;
    int64_t out[MAX_VTX];
    int64_t wgt[MAX_VTX];

    stinger_gather_successors(graph, from, &outlen, out, wgt, NULL, NULL, NULL, MAX_VTX);
    fail_unless(outlen >= 0, "Length of output array %ld < 0", outlen);

    edgemap_t graph_edges_copy;
    edgemap_t::iterator iter;
    for (iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      graph_edges_copy.insert( pair<edge_t,edge_info_t>((*iter).first, (*iter).second) );
    }

    for (edgemap_t::iterator iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      edge_t edge = (*iter).first;
      if (edge.first != from || (edge.first == edge.second) || (*iter).second.type != 0) {
	graph_edges_copy.erase(edge);
      }
    }

    for (int j = 0; j < outlen; ++j) {

      edgemap_t::iterator iter = graph_edges_copy.find(edge_t(from, out[j]));
      fail_unless(iter != graph_edges_copy.end(), "Non Existent edge returned");

      edge_t edge = (*iter).first;
      edge_info_t einfo = (*iter).second;

      fail_unless(edge.first == from, "Internal Error");
      fail_unless(edge.second == out[j], "Internal Error");

      fail_unless(einfo.wgt == wgt[j], "Weight Does not match, Expected %ld, Returned %ld", einfo.wgt, wgt[j]);

      graph_edges_copy.erase((*iter).first);
    }

    fail_unless(graph_edges_copy.size() == 0, "Edges Remaining");
  }

}
END_TEST

START_TEST(test_stinger_gather_typed_successors_serial)
{
  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  for (int i = 0; i < NV; ++i) {
    if (i != 0) continue;
    int from = i;
    size_t outlen;
    int64_t out[MAX_VTX];

    stinger_gather_typed_successors_serial(graph, 0, from, &outlen, out, MAX_VTX);
    fail_unless(outlen >= 0, "Length of output array %ld < 0", outlen);

    edgemap_t graph_edges_copy;
    edgemap_t::iterator iter;
    for (iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      graph_edges_copy.insert( pair<edge_t,edge_info_t>((*iter).first, (*iter).second) );
    }

    for (edgemap_t::iterator iter = graph_edges->begin(); iter != graph_edges->end(); ++iter) {
      edge_t edge = (*iter).first;
      if (edge.first != from || (*iter).second.type != 0) {
	graph_edges_copy.erase(edge);
      }
    }

    for (int j = 0; j < outlen; ++j) {

      edgemap_t::iterator iter = graph_edges_copy.find(edge_t(from, out[j]));
      fail_unless(iter != graph_edges_copy.end(), "Non Existent edge returned");

      edge_t edge = (*iter).first;

      fail_unless(edge.first == from, "Internal Error");
      fail_unless(edge.second == out[j], "Internal Error");

      graph_edges_copy.erase((*iter).first);
    }

    fail_unless(graph_edges_copy.size() == 0, "Edges Remaining");
  }

}
END_TEST

START_TEST(test_stinger_to_sorted_csr)
{
  for (int i = 0; i < NE; ++i) {
    stinger_insert_edge(graph, graph_types[i], 
	graph_source[i], graph_dest[i],
	graph_wgts[i], graph_timestamps[i]);
  }

  int64_t *off;
  int64_t *ind;
  int64_t *wgt;
  stinger_to_sorted_csr(graph, NV, &off, &ind, &wgt, NULL, NULL, NULL);
  
  struct stinger* graph_copy = graph;
  graph = stinger_new();

  for (int64_t v = 0; v < NV; ++v) {
    int64_t prev = -1;
    for (int64_t offset = off[v]; offset < off[v+1]; ++offset) {
      fail_unless(ind[offset] >= prev, "Ordering Not Correct for source vertex %ld: %ld < %ld", v, ind[offset], prev);
      prev = ind[offset];
    }
  }

  stinger_set_initial_edges(graph, NV, 0, off, ind, wgt, NULL, NULL, 0);
  check_edges_exist(CHECK_ALL & ~CHECK_TIMESTAMPS);
  //stinger_free(graph);
  graph = graph_copy;
}
END_TEST

  Suite *
stinger_basic_suite (void)
{
  Suite *s = suite_create ("Stinger-Basic");

  for (int i = 0; i < graph_num_tests_untyped; ++i) {

    TCase *tc_graph = tcase_create (graph_name[i]);

    tcase_add_checked_fixture(tc_graph, graph_setup[i], graph_cleanup[i]);

    int num_tests = test_names.size() - 1;
    for (int j = 1; j < num_tests; ++j) {
      if (tests_to_run[j]) {
	map<TEST_ID, string>::iterator iter = test_ids.find((TEST_ID)j);
	if (iter != test_ids.end())
	  _tcase_add_test(tc_graph, test_funcs[j], (*iter).second.c_str(), 0, 0, 0, 1);
      }
    }
    suite_add_tcase(s, tc_graph);
  }

  for (int i = 0; i < num_random_tests; ++i) {

    TCase *tc_graph = tcase_create ("Random Graph Test");

    tcase_add_checked_fixture(tc_graph, graph_random_setup, graph_random_cleanup);
    int num_tests = test_names.size() - 1;
    for (int j = 1; j < num_tests; ++j) {
      if (tests_to_run[j]) {
	map<TEST_ID, string>::iterator iter = test_ids.find((TEST_ID)j);
	if (iter != test_ids.end())
	  _tcase_add_test(tc_graph, test_funcs[j], (*iter).second.c_str(), 0, 0, 0, 1);
      }
    }

    suite_add_tcase(s, tc_graph);
  }

  return s;
}

int main(int argc, char** argv)
{
  init_test_registry();
  int opt;
  while ((opt = getopt(argc, argv, "at:")) != -1) {
    switch(opt) {
      case 'a': // All tests
	memset(tests_to_run, 0xFF, 1000*sizeof(char));
	break;
      case 't':
	char* tname = optarg;
	TEST_ID id = find_test_id(tname);
	tests_to_run[id] = 0xFF;
	break;
    }   
  }
  //  srand(time(0));
  int number_failed;
  Suite *s = stinger_basic_suite();
  SRunner *sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

void init_test_registry()
{
  REGISTER(INVALID, NULL);
  REGISTER(STINGER_INSERT_EDGE, stinger_insert_edge);
  REGISTER(STINGER_INSERT_EDGE_PAIR, stinger_insert_edge_pair);
  REGISTER(STINGER_REMOVE_EDGE, stinger_remove_edge);
  REGISTER(STINGER_INCR_EDGE_PAIR, stinger_incr_edge_pair);
  REGISTER(STINGER_REMOVE_EDGE_PAIR, stinger_remove_edge_pair);
  REGISTER(STINGER_INCR_EDGE, stinger_incr_edge);
  REGISTER(STINGER_TOTAL_EDGES, stinger_total_edges);
  REGISTER(STINGER_REMOVE_AND_INSERT_EDGES, stinger_remove_and_insert_edges);
  REGISTER(STINGER_REMOVE_AND_INSERT_EDGES_BATCH, stinger_remove_and_insert_edges_batch);
  REGISTER(STINGER_GATHER_TYPED_SUCCESSORS, stinger_gather_typed_successors);
  REGISTER(STINGER_GATHER_SUCCESSORS, stinger_gather_successors);
  REGISTER(STINGER_GATHER_TYPED_SUCCESSORS_SERIAL, stinger_gather_typed_successors_serial);
  REGISTER(STINGER_TO_SORTED_CSR, stinger_to_sorted_csr);
}

TEST_ID find_test_id(char* test_name) {

  assert(test_name != NULL);
  string str = string(test_name);
  std::transform(str.begin(), str.end(),str.begin(), ::toupper);

  map<string, TEST_ID>::iterator test_iter = test_names.find(str);
  if (test_iter == test_names.end()) {
    return TEST_INVALID;
  }

  return (*test_iter).second;
}

void test_NULL(int a)
{
}

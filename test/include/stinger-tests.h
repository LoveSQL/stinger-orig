#ifndef STINGER_TESTS_H
#define STINGER_TESTS_H


extern "C"
{
	#include <stinger.h>
	#include <stinger-utils.h>
}

#include <map>

#define MAX_VTX 100


struct edge_info_t
{
	int64_t type;
	int64_t wgt;
	int64_t timestamp;

	edge_info_t(int64_t t, int64_t w, int64_t ts) {
		type = t;
		wgt = w;
		timestamp = ts;
	}

	edge_info_t& operator=(const edge_info_t& ei) {
		type = ei.type;
		wgt = ei.wgt;
		timestamp = ei.timestamp;
		return (*this);
	}
};

#define PRINT_ARRAY(arr,beg,end,...) do\
{\
	printf(__VA_ARGS__);\
	int64_t __i_pa;\
	for (__i_pa = beg; __i_pa <= end; ++__i_pa){\
		printf("%ld ", arr[__i_pa]);\
	}\
	printf("\n");\
} while(0)

#define PRINT_STINGER_EDGES(__STINGER, __TYPE, ...) do {\
	STINGER_FORALL_EDGES_BEGIN(__STINGER, __TYPE)\
		printf(__VA_ARGS__);\
		printf("%ld %ld", STINGER_EDGE_SOURCE, STINGER_EDGE_DEST);\
	STINGER_FORALL_EDGES_END();\
} while(0)

#endif

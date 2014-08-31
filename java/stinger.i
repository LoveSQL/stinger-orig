%module stinger
%{
#include "stdint.h"
#include "stinger.h"
#include "stinger-utils.h"
#include "stinger-physmap.h"
#include "stinger-atomics.h"
#include "stinger-iterator.h"
%}
#define XMTI    
#define const   

%include "stdint.i"

%include "stinger.h"
%include "stinger-utils.h"
%include "stinger-physmap.h"
%include "stinger-atomics.h"
%include "stinger-iterator.h"
%include "carrays.i"
%array_class(uint64_t, uint64Array);
%array_class(int64_t, int64Array);

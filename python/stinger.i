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

/* TRICKY STDINT STUFF */
%include <swigarch.i>

/* Exact integral types.  */

/* Signed.  */

#define int8_t signed char  
#define int16_t short int  
#define int32_t int    
#if defined(SWIGWORDSIZE64)
#define int64_t long int   
#else
#define int64_t long long int  
#endif

/* Unsigned.  */
#define uint8_t unsigned char  
#define uint16_t unsigned short int 
#define uint32_t unsigned int   
#if defined(SWIGWORDSIZE64)
#define uint64_t unsigned long int
#else
#define uint64_t unsigned long long int 
#endif


/* Small types.  */

/* Signed.  */
#define int_least8_t signed char  
#define int_least16_t short int  
#define int_least32_t int    
#if defined(SWIGWORDSIZE64)
#define int_least64_t long int   
#else
#define int_least64_t long long int  
#endif

/* Unsigned.  */
#define uint_least8_t unsigned char  
#define uint_least16_t unsigned short int 
#define uint_least32_t unsigned int   
#if defined(SWIGWORDSIZE64)
#define uint_least64_t unsigned long int
#else
#define uint_least64_t unsigned long long int 
#endif


/* Fast types.  */

/* Signed.  */
#define int_fast8_t signed char  
#if defined(SWIGWORDSIZE64)
#define int_fast16_t long int   
#define int_fast32_t long int   
#define int_fast64_t long int   
#else
#define int_fast16_t int    
#define int_fast32_t int    
#define int_fast64_t long long int  
#endif

/* Unsigned.  */
#define uint_fast8_t unsigned char  
#if defined(SWIGWORDSIZE64)
#define uint_fast16_t unsigned long int
#define uint_fast32_t unsigned long int
#define uint_fast64_t unsigned long int
#else
#define uint_fast16_t unsigned int   
#define uint_fast32_t unsigned int   
#define uint_fast64_t unsigned long long int 
#endif


/* Types for `void *' pointers.  */
#if defined(SWIGWORDSIZE64)
#define intptr_t long int   
#define uintptr_t unsigned long int
#else
#define intptr_t int    
#define uintptr_t unsigned int   
#endif


/* Largest integral types.  */
#if defined(SWIGWORDSIZE64)
#define intmax_t long int   
#define uintmax_t unsigned long int
#else
#define intmax_t long long int  
#define uintmax_t unsigned long long int 
#endif

/* END TRICKY STDINT STUFF */


%include "stinger.h"
%include "stinger-utils.h"
%include "stinger-physmap.h"
%include "stinger-atomics.h"
%include "carrays.i"
%include "stinger-iterator.h"
%array_class(uint64_t, uint64Array);
%array_class(int64_t, int64Array);

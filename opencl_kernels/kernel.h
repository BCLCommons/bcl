#ifndef KERNEL_H_
#define KERNEL_H_

#define PRECISION float
#define __global
#define __local
#define __kernel

#define true  1
#define false 0

#define uint   unsigned int
#define size_t uint

#define get_local_size( x)  1
#define get_global_size( x) 1

#define get_local_id( x)    1
#define get_global_id( x)   1
#define get_group_id( x)    1

#define get_num_groups( x)  1

#define barrier( x)

#define float2              union{ float x; float y;}
#define float4              union{ float x; float y; float z; float w;}

#define INFINITY            999999.99

#define min( x, y)          x
#define max( x, y)          x

#define isequal( x, y)        ( x == y)
#define isnotequal( x, y)     ( x != y)
#define isgreater( x, y)      ( x >  y)
#define isgreaterequal( x, y) ( x >= y)
#define isless( x, y)         ( x <  y)
#define islessequal( x, y)    ( x <= y)
#define isfinite( x)          false
#define isinf( x)             false
#define isnan( x)             false
#define isnormal( x)          false
#define isordered( x)         false
#define isunordered( x)       false
#define signbit( x)           false
#define any( x)               false
#define all( x)               false
#define bitselect( x)         false
#define select( x)            false

#define native_sqrt( x)      sqrt( x)
#define native_recip( x)     1.0 / x
#define native_exp( x)       exp( x)
#define native_divide( x, y) x / y

#define dot( x)             x
#define mad( x, y, z)       x

#define pown( x, n)         x

// macros that should be defined
#define FLT_DIG 6
#define FLT_MANT_DIG  24
#define FLT_MAX_10_EXP  +38
#define FLT_MAX_EXP +128
#define FLT_MIN_10_EXP  -37
#define FLT_MIN_EXP -125
#define FLT_RADIX 2
#define FLT_MAX 0x1.fffffep127f
#define FLT_MIN 0x1.0p-126f
#define FLT_EPSILON 0x1.0p-23f
#define CHAR_BIT  8
#define CHAR_MAX  SCHAR_MAX
#define CHAR_MIN  SCHAR_MIN
#define INT_MAX 2147483647
#define INT_MIN (-2147483647- 1)
#define LONG_MAX  0x7fffffffffffffffL
#define LONG_MIN  (-0x7fffffffffffffffL- 1)
#define SCHAR_MAX 127
#define SCHAR_MIN (-127 - 1)
#define SHRT_MAX  32767
#define SHRT_MIN  (-32767- 1)
#define UCHAR_MAX 255
#define USHRT_MAX 65535
#define UINT_MAX  0xffffffff
#define ULONG_MAX 0xffffffffffffffffUL

// double
#define DBL_DIG 15
#define DBL_MANT_DIG  53
#define DBL_MAX_10_EXP  +308
#define DBL_MAX_EXP +1024
#define DBL_MIN_10_EXP  -307
#define DBL_MIN_EXP -1021
#define DBL_MAX 0x1.fffffffffffffp1023
#define DBL_MIN 0x1.0p-1022
#define DBL_EPSILON 0x1.0p-52

#endif /* KERNEL_H_ */

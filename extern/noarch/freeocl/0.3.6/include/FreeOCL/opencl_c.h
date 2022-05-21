/*
	FreeOCL - a free OpenCL implementation for CPU
	Copyright (C) 2011  Roland Brochard

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>
*/
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_H__

#include "FreeOCL/config.h"
#include <cstddef>

#define __OPENCL_VERSION__	120
#define CL_VERSION_1_0		100
#define CL_VERSION_1_1		110
#define CL_VERSION_1_2		120

#define MAXFLOAT		0x1.fffffep127f
//#define HUGE_VALF
//#define INFINITY
//#define NAN

#define FLT_DIG			6
#define FLT_MANT_DIG	24
#define FLT_MAX_10_EXP	+38
#define FLT_MAX_EXP		+128
#define FLT_MIN_10_EXP	-37
#define FLT_MIN_EXP		-125
#define FLT_RADIX		2
#define FLT_MAX			0x1.fffffep127f
#define FLT_MIN			0x1.0p-126f
#define FLT_EPSILON		0x1.0p-23f

#define M_E_F			2.7182818284590452354f	/* e */
#define M_LOG2E_F		1.4426950408889634074f	/* log_2 e */
#define M_LOG10E_F		0.43429448190325182765f	/* log_10 e */
#define M_LN2_F			0.69314718055994530942f	/* log_e 2 */
#define M_LN10_F		2.30258509299404568402f	/* log_e 10 */
#define M_PI_F			3.14159265358979323846f	/* pi */
#define M_PI_2_F		1.57079632679489661923f	/* pi/2 */
#define M_PI_4_F		0.78539816339744830962f	/* pi/4 */
#define M_1_PI_F		0.31830988618379067154f	/* 1/pi */
#define M_2_PI_F		0.63661977236758134308f	/* 2/pi */
#define M_2_SQRTPI_F	1.12837916709551257390f	/* 2/sqrt(pi) */
#define M_SQRT2_F		1.41421356237309504880f	/* sqrt(2) */
#define M_SQRT1_2_F		0.70710678118654752440f	/* 1/sqrt(2) */

#define CHAR_BIT	8
#define CHAR_MAX	SCHAR_MAX
#define CHAR_MIN	SCHAR_MIN
#define INT_MAX		2147483647
#define INT_MIN		(-2147483647 - 1)
#define LONG_MAX	0x7fffffffffffffffL
#define LONG_MIN	(-0x7fffffffffffffffL - 1)
#define SCHAR_MAX	127
#define SCHAR_MIN	(-127 - 1)
#define SHRT_MAX	32767
#define SHRT_MIN	(-32767 - 1)
#define UCHAR_MAX	255
#define USHRT_MAX	65535
#define UINT_MAX	0xffffffff
#define ULONG_MAX	0xffffffffffffffffUL

#define __kernel
#define kernel	__kernel

#define FLOAT2	__m64
#define FLOAT4	__m128
#define packed	__packed__

#ifdef __size_t
#undef __size_t
#endif

// Built-in scalar types
#ifdef __GNUG__

#ifdef __i386__
#include <stdint.h>
typedef char			__char;
typedef uint8_t         __uchar;
typedef int16_t         __short    __attribute__((aligned(2)));
typedef uint16_t        __ushort   __attribute__((aligned(2)));
typedef int32_t         __int      __attribute__((aligned(4)));
typedef uint32_t        __uint     __attribute__((aligned(4)));
typedef int64_t         __long     __attribute__((aligned(8)));
typedef uint64_t        __ulong    __attribute__((aligned(8)));

typedef float           __float    __attribute__((aligned(4)));
typedef double          __double   __attribute__((aligned(8)));

typedef bool			__bool;
typedef __ulong			__size_t;

#elif defined __x86_64__
#include <stdint.h>
typedef char			__char;
typedef uint8_t         __uchar;
typedef int16_t         __short    __attribute__((aligned(2)));
typedef uint16_t        __ushort   __attribute__((aligned(2)));
typedef int32_t         __int      __attribute__((aligned(4)));
typedef uint32_t        __uint     __attribute__((aligned(4)));
typedef int64_t         __long     __attribute__((aligned(8)));
typedef uint64_t        __ulong    __attribute__((aligned(8)));

typedef float           __float    __attribute__((aligned(4)));
typedef double          __double   __attribute__((aligned(8)));

typedef bool			__bool;
typedef __ulong			__size_t;

#else
typedef bool			__bool;
typedef char			__char;
typedef unsigned char	__uchar;
typedef short			__short;
typedef unsigned short	__ushort;
typedef int				__int;
typedef unsigned int	__uint;
typedef long			__long;
typedef unsigned long	__ulong;
typedef float			__float;
typedef double			__double;
typedef __ulong			__size_t;
#endif

#else
typedef bool			__bool;
typedef char			__char;
typedef unsigned char	__uchar;
typedef short			__short;
typedef unsigned short	__ushort;
typedef int				__int;
typedef unsigned int	__uint;
typedef long			__long;
typedef unsigned long	__ulong;
typedef float			__float;
typedef double			__double;
typedef __ulong			__size_t;
#endif

struct half
{
	__ushort v;

	inline operator float() const
	{
		union
		{
			__uint i;
			float f;
		} u;
		const __uint s = (v & 0x8000U) << 16;
		const __uint e = (v & 0x7C00U) >> 10;
		const __uint m = v & 0x03FFU;
		switch(e)
		{
		case 0:
			u.i = s;	// flush denormalized numbers to 0 for now
			break;
		case 0x1F:
			u.i = s				// signbit
				  | 0x7F800000U				// exponent
				  | (m << 13);	// significand
			break;
		default:
			u.i = s				// signbit
				  | ((e + 126 - 15) << 23)	// exponent
				  | (m << 13);	// significand
		}
		return u.f;
	}

	static inline half from_float(float f)
	{
		union
		{
			__uint i;
			float f;
		} u;
		u.f = f;

		half h;
		const __uint s = (u.i & 0x80000000U) >> 16;
		const __uint e = (u.i & 0x7F800000U) >> 23;
		const __uint m = (u.i & 0x007FFFFFU) >> 13;
		switch(e)
		{
		case 0:
			h.v = s;
			break;
		case 0xFF:
			h.v = s | 0x7C00U | m;
			break;
		default:
			h.v = s | (e << 10) | m;
		}
		return h;
	}
};
#ifdef __GNUG__
typedef half __half     __attribute__((aligned(2)));
#else
typedef half __half;
#endif

template<class A, class B>	struct __right {	typedef B	type;	};
template<class A, class B>	struct __match {	enum { value = false };	};
template<class A>	struct __match<A, A> {	enum { value = true };	};

// A small template to filter scalar types allowed for arithmetic operations
template<typename S>	struct __scalar;
template<>	struct __scalar<__char>		{	typedef __char	type;	};
template<>	struct __scalar<__short>	{	typedef __short	type;	};
template<>	struct __scalar<__int>		{	typedef __int	type;	};
template<>	struct __scalar<__long>		{	typedef __long	type;	};
template<>	struct __scalar<__uchar>	{	typedef __uchar	type;	};
template<>	struct __scalar<__ushort>	{	typedef __ushort	type;	};
template<>	struct __scalar<__uint>		{	typedef __uint	type;	};
template<>	struct __scalar<__ulong>	{	typedef __ulong	type;	};
template<>	struct __scalar<__float>	{	typedef __float	type;	};
template<>	struct __scalar<__double>	{	typedef __double	type;	};

template<typename B>	struct __sint_type_of_same_size;
template<>	struct __sint_type_of_same_size<__char>		{	typedef __char	type;	};
template<>	struct __sint_type_of_same_size<__uchar>	{	typedef __char	type;	};
template<>	struct __sint_type_of_same_size<__short>	{	typedef __short	type;	};
template<>	struct __sint_type_of_same_size<__ushort>	{	typedef __short	type;	};
template<>	struct __sint_type_of_same_size<__int>		{	typedef __int	type;	};
template<>	struct __sint_type_of_same_size<__uint>		{	typedef __int	type;	};
template<>	struct __sint_type_of_same_size<__long>		{	typedef __long	type;	};
template<>	struct __sint_type_of_same_size<__ulong>	{	typedef __ulong	type;	};
template<>	struct __sint_type_of_same_size<__float>	{	typedef __int	type;	};
template<>	struct __sint_type_of_same_size<__double>	{	typedef __long	type;	};

template<bool, class T> struct __if;
template<class T> struct __if<true, T>	{	typedef T type;	};

template<typename T, int Id>
static inline T &__create_local()
{
	static __thread T v;
	return v;
}
template<typename T, int Id>
static inline T &__create_local(const T &init)
{
	static __thread T v;
	memcpy(&v, &init, sizeof(T));
	return v;
}

#include "vectors.h"

template<typename S, typename V = void>	struct __igentype;
template<>	struct __igentype<__char>	{	typedef __char	type;	};
template<>	struct __igentype<__short>	{	typedef __short	type;	};
template<>	struct __igentype<__int>	{	typedef __int	type;	};
template<>	struct __igentype<__long>	{	typedef __long	type;	};
template<>	struct __igentype<__uchar>	{	typedef __uchar	type;	};
template<>	struct __igentype<__ushort>	{	typedef __ushort	type;	};
template<>	struct __igentype<__uint>	{	typedef __uint	type;	};
template<>	struct __igentype<__ulong>	{	typedef __ulong	type;	};
template<class V>
struct __igentype<V, typename __right<__igentype<typename __vector<V>::base_type>, void>::type >
{
	typedef V	type;
};

#include "simdopts.h"

#include "workitem.h"
#include "math.h"
#include "integer.h"
#include "common.h"
#include "geometric.h"
#include "relational.h"
#include "vloadstore.h"
#include "sync.h"
#include "memfence.h"
#include "asyncmemop.h"
#include "atomic.h"
#include "vmisc.h"
#include "converters.h"
#include "imgreadwrite.h"
#include "printf.h"

#undef FLOAT2
#undef FLOAT4

#define __global
#define __local
#define __constant	const
#define __private
#define global		__global
#define local		__local
#define constant	__constant
#define private		__private

#ifdef __GNUG__
#define restrict	__restrict__
#else
#define restrict
#endif

#define __read_only	const
#define __write_only
#define __read_write
#define read_only	__read_only
#define write_only	__write_only
#define read_write	__read_write

#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1_I
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1_UI
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2_I
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2_UI
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2I
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PI
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PB
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2S
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2S_
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3_I
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3_UI
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3SS
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3_S
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3S_
#undef VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3SS_
#undef CHECK

// debugging helpers

template<int N, typename T, typename I>
static inline T &__check_lookup_bounds(T *ptr, const I &idx, const char *ref)
{
	if (idx < 0 || idx >= N)
		printf("warning: (%d,%d,%d|%d,%d,%d) %d is out of array range [0 - %d[ in %s\n",
			   get_global_id(0), get_global_id(1), get_global_id(2),
			   get_local_id(0), get_local_id(1), get_local_id(2),
			   idx, N, ref);
	return ptr[idx];
}

template<int N, typename T, typename I>
static inline const T &__check_lookup_bounds(const T *ptr, const I &idx, const char *ref)
{
	if (idx < 0 || idx >= N)
		printf("warning: (%d,%d,%d|%d,%d,%d) %d is out of array range [0 - %d[ in %s\n",
			   get_global_id(0), get_global_id(1), get_global_id(2),
			   get_local_id(0), get_local_id(1), get_local_id(2),
			   idx, N, ref);
	return ptr[idx];
}

#endif

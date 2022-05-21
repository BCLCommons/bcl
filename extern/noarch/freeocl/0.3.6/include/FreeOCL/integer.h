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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_INTEGER_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_INTEGER_H__
#include "vectors.h"

// Built-in integer functions

static inline __long clamp(__long, __long, __long);
static inline __ulong clamp(__ulong, __ulong, __ulong);

// for __int
static inline __uint __abs(const __int &x)	{	return x >= 0 ? x : -x;	}
#define abs(x)	__abs(x)
static inline __uint abs_diff(__int x, __int y)	{	return x > y ? x - y : y - x;	}
static inline __int add_sat(__int x, __int y)
{	return __int(clamp(__long(x) + __long(y), -0x80000000, 0x7FFFFFFF));	}

static inline __int hadd(__int x, __int y)
{	return __int((__long(x) + __long(y)) >> 1);	}

static inline __int rhadd(__int x, __int y)
{	return __int((__long(x) + __long(y) + 1) >> 1);	}

static inline __int max(__int, __int);
static inline __int min(__int, __int);
static inline __int clamp(__int x, __int minval, __int maxval)	{	return min(max(x, minval), maxval);	}
static inline __int clz(__int x)
{
	const __uint y = (__uint)x;
	__int c = ((y >> 16) == 0) << 4;
	c += ((y >> (24 - c)) == 0) << 3;
	c += ((y >> (28 - c)) == 0) << 2;
	c += ((y >> (30 - c)) == 0) << 1;
	c += ((y >> (31 - c)) == 0);
	return c;
}
static inline __int mul_hi(__int, __int);
static inline __int mad_hi(__int a, __int b, __int c)	{	return mul_hi(a, b) + c;	}
static inline __int mad_sat(__int a, __int b, __int c)
{	return __int(clamp(__long(a) * __long(b) + __long(c), (__long)-0x80000000, (__long)0x7FFFFFFF));	}

static inline __int max(__int x, __int y)	{	return x > y ? x : y;	}
static inline __uint max(__uint x, __uint y)	{	return x > y ? x : y;	}
static inline __int min(__int x, __int y)	{	return x < y ? x : y;	}
static inline __uint min(__uint x, __uint y)	{	return x < y ? x : y;	}
static inline __int mul_hi(__int x, __int y)	{	return __int((long long)(x) * (long long)(y) >> 32);	}
static inline __int rotate(__int v, __int i)	{	return (v << i) | (v >> (32 - i));	}

static inline __int sub_sat(__int x, __int y)
{	return __int(clamp(__long(x) - __long(y), (__long)-0x80000000, (__long)0x7FFFFFFF));	}

// for __long
static inline __long max(__long, __long);
static inline __long min(__long, __long);
static inline __long clamp(__long x, __long minval, __long maxval)	{	return min(max(x, minval), maxval);	}

static inline __long max(__long x, __long y)	{	return x > y ? x : y;	}
static inline __ulong max(__ulong x, __ulong y)	{	return x > y ? x : y;	}
static inline __long min(__long x, __long y)	{	return x < y ? x : y;	}
static inline __ulong min(__ulong x, __ulong y)	{	return x < y ? x : y;	}
static inline __ulong clamp(__ulong x, __ulong minval, __ulong maxval)	{	return min(max(x, minval), maxval);	}

// upsample functions
static inline __short upsample(__char hi, __uchar lo)	{	return ((__short)hi << 8) | lo;	}
static inline __ushort upsample(__uchar hi, __uchar lo)	{	return ((__ushort)hi << 8) | lo;	}
static inline __int upsample(__short hi, __ushort lo)	{	return ((__int)hi << 16) | lo;	}
static inline __uint upsample(__ushort hi, __ushort lo)	{	return ((__uint)hi << 16) | lo;	}
static inline __long upsample(__int hi, __uint lo)	{	return ((long long)hi << 32) | lo;	}
static inline __ulong upsample(__uint hi, __uint lo)	{	return ((long long)hi << 32) | lo;	}

static inline __ulong popcount(__ulong x)
{
#ifdef __GNUG__
	return __builtin_popcount(x);
#else
	const __ulong m1 = 0x5555555555555555;
	const __ulong m2 = 0x3333333333333333;
	const __ulong m4 = 0x0f0f0f0f0f0f0f0f;
	x -= (x >> 1) & m1;
	x = (x & m2) + ((x >> 2) & m2);
	x = (x + (x >> 4)) & m4;
	return (x * h01) >> 56;
#endif
}

// fast integer built-in functions
static inline __int mad24(__int x, __int y, __int z)	{	return x * y + z;	}
static inline __int mul24(__int x, __int y)	{	return x * y;	}

VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3_I(clamp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1_I(abs)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1_I(popcount)

#endif

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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_RELATIONAL_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_RELATIONAL_H__
// undefining C99 macros for some of these functions (they are functions in OpenCL C, not macros)
template<class F>
static inline __int __isfinite(F x)	{	return isfinite(x);	}
template<class F>
static inline __int __isinf(F x)	{	return isinf(x);	}
template<class F>
static inline __int __isnan(F x)	{	return isnan(x);	}
template<class F>
static inline __int __isnormal(F x)	{	return isnormal(x);	}
template<class F>
static inline __int __isunordered(F x, F y)	{	return isunordered(x, y);	}
template<class F>
static inline __int __signbit(F x)	{	return signbit(x);	}

#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef isunordered
#undef signbit

// Built-in relational functions
// for scalar types
template<typename S>
static inline __int isequal(S x, typename __scalar<S>::type y)	{	return x == y;	}
template<typename S>
static inline __int isnotequal(S x, typename __scalar<S>::type y)	{	return x != y;	}
template<typename S>
static inline __int isgreater(S x, typename __scalar<S>::type y)	{	return x > y;	}
template<typename S>
static inline __int isgreaterequal(S x, typename __scalar<S>::type y)	{	return x >= y;	}
template<typename S>
static inline __int isless(S x, typename __scalar<S>::type y)	{	return x < y;	}
template<typename S>
static inline __int islessequal(S x, typename __scalar<S>::type y)	{	return x <= y;	}
template<typename S>
static inline __int islessgreater(S x, typename __scalar<S>::type y)	{	return x < y || x > y;	}

static inline __int isfinite(__float x)	{	return __isfinite(x);	}
static inline __int isinf(__float x)	{	return __isinf(x);	}
static inline __int isnan(__float x)	{	return __isnan(x);	}
static inline __int isnormal(__float x)	{	return __isnormal(x);	}

template<typename S>
static inline __int isordered(S x, typename __scalar<S>::type y)	{	return !__isordered(x, y);	}

template<typename S>
static inline __int isunordered(S x, typename __scalar<S>::type y)	{	return __isunordered(x, y);	}

static inline __int signbit(__float x)	{	return __signbit(x);	}

// both scalars and vectors
template<typename S>
static inline S bitselect(S a, S b, S c)
{
	S r;
	for(__size_t i = 0 ; i < sizeof(S) ; ++i)
		(&reinterpret_cast<char&>(r))[i]
				= ((&reinterpret_cast<char&>(b))[i] & (&reinterpret_cast<char&>(c))[i])
		| ((&reinterpret_cast<char&>(a))[i] & ~(&reinterpret_cast<char&>(c))[i]);
	return r;
}

template<typename S, typename I>
static inline typename __scalar<S>::type select(S a, S b, I c)
{	return c ? b : a;	}

// for vector types

template<typename I>
static inline __int any(I x)
{
	const __size_t mask = 1U << (sizeof(typename __vector<I>::base_type) * 8 - 1);
	for(__size_t i = 0 ; i < __vector<I>::components ; ++i)
		if (x.v[i] & mask)
			return 1;
	return 0;
}

template<typename I>
static inline __int all(I x)
{
	const __size_t mask = 1U << (sizeof(typename __vector<I>::base_type) * 8 - 1);
	for(__size_t i = 0 ; i < __vector<I>::components ; ++i)
		if (!(x.v[i] & mask))
			return 0;
	return 1;
}

template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isequal(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] == y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isnotequal(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] != y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isgreater(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] > y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isgreaterequal(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] >= y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isless(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] < y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		islessequal(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = x.v[i] <= y.v[i] ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		islessgreater(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = (x.v[i] < y.v[i] || x.v[i] > y.v[i]) ? -1 : 0;
	return ret;
}

template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isfinite(S x)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isfinite(x.v[i]) ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isinf(S x)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isinf(x.v[i]) ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isnan(S x)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isnan(x.v[i]) ? -1 : 0;
	return ret;
}
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isnormal(S x)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isnormal(x.v[i]) ? -1 : 0;
	return ret;
}

template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isordered(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isunordered(x.v[i], y.v[i]) ? 0 : -1;
	return ret;
}

template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		isunordered(S x, typename __vector<S>::type y)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __isunordered(x.v[i], y.v[i]) ? -1 : 0;
	return ret;
}

// This is implemented in math.h
template<typename S>
static inline typename __vector_type<__int, __vector<S>::components>::type
		signbit(S x)
{
	typename __vector_type<__int, __vector<S>::components>::type ret;
	for(__size_t i = 0 ; i < __vector<S>::components ; ++i)
		ret.v[i] = __signbit(x.v[i]) ? -1 : 0;
	return ret;
}

template<typename V, typename I>
static inline typename __vector<V>::type select(V a, V b, I c)
{
	V ret;
	const __size_t mask = 1U << (sizeof(typename __vector<I>::base_type) * 8 - 1);
	for(__size_t i = 0 ; i < __vector<V>::components ; ++i)
		ret.v[i] = (c.v[i] & mask) ? b.v[i] : a.v[i];
	return ret;
}

#endif

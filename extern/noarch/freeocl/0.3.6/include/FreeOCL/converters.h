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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_CONVERTERS_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_CONVERTERS_H__

template<class T> struct __type_range;
template<> struct __type_range<__char>		{	enum { __min = -128, __max = 127 };	};
template<> struct __type_range<__uchar>		{	enum { __min = 0, __max = 255 };	};
template<> struct __type_range<__short>		{	enum { __min = -32768, __max = 32767 };	};
template<> struct __type_range<__ushort>	{	enum { __min = 0, __max = 65536 };	};
template<> struct __type_range<__int>		{	enum { __min = -0x10000000, __max = 0x7FFFFFFF };	};
template<> struct __type_range<__uint>		{	enum { __min = 0, __max = 0xFFFFFFFFU };	};
template<> struct __type_range<__long>		{	enum { __min = -0x1000000000000000L, __max = 0x7FFFFFFFFFFFFFFFL };	};
template<> struct __type_range<__ulong>		{	enum { __min = 0, __max = 0xFFFFFFFFFFFFFFFFLU };	};

// Built-in type conversion functions
#define CONVERTER(X)\
template<class Scalar>\
static inline typename __if< !__match<typename __scalar<Scalar>::type, __##X>::value, __##X>::type convert_##X(const Scalar &v)	{	return (__##X)v;	}\
template<class Scalar>\
static inline const typename __if< __match<Scalar, __##X>::value, __##X>::type &convert_##X(const Scalar &v)	{	return v;	}

#define CONVERTER_SAT(X)\
template<class Scalar>\
static inline typename __if< !__match<typename __scalar<Scalar>::type, __##X>::value, __##X>::type convert_##X##_sat(const Scalar &v)\
{\
	return (__##X)(v > __type_range<__##X>::__max\
			   ? __type_range<__##X>::__max\
			   : v < __type_range<__##X>::__min\
				 ? __type_range<__##X>::__min\
				 : v);\
}\
template<class Scalar>\
static inline const typename __if< __match<Scalar, __##X>::value , __##X>::type &convert_##X##_sat(const Scalar &v)	{	return v;	}

CONVERTER(char)
CONVERTER(uchar)
CONVERTER(short)
CONVERTER(ushort)
CONVERTER(int)
CONVERTER(uint)
CONVERTER(long)
CONVERTER(ulong)
CONVERTER(float)
CONVERTER(double)

CONVERTER_SAT(char)
CONVERTER_SAT(uchar)
CONVERTER_SAT(short)
CONVERTER_SAT(ushort)
CONVERTER_SAT(int)
CONVERTER_SAT(uint)
CONVERTER_SAT(long)
CONVERTER_SAT(ulong)

#undef CONVERTER
#undef CONVERTER_SAT

#define CONVERTER(X, N)\
template<class V>\
	static inline typename __if<__vector<V>::components == N && !__match<__##X##N, V>::value, __##X##N>::type convert_##X##N(const V &v)\
{\
	__##X##N ret;\
	for(size_t i = 0 ; i < N ; ++i)\
		ret.v[i] = (__##X)v.v[i];\
	return ret;\
}\
template<class V>\
	static inline const typename __if<__match<__##X##N, V>::value, __##X##N>::type &convert_##X##N(const V &v)	{	return v;	}
#define CONVERTER_SAT(X, N)\
template<class V>\
	static inline typename __if<__vector<V>::components == N && !__match<__##X##N, V>::value, __##X##N>::type convert_##X##N##_sat(const V &v)\
{\
	__##X##N ret;\
	for(size_t i = 0 ; i < N ; ++i)\
		ret.v[i] = convert_##X##_sat(v.v[i]);\
	return ret;\
}\
template<class V>\
	static inline const typename __if<__match<__##X##N, V>::value, __##X##N>::type &convert_##X##N##_sat(const V &v)	{	return v;	}

#define CONVERTER_VECTORS(X)\
	CONVERTER(X, 2)\
	CONVERTER(X, 3)\
	CONVERTER(X, 4)\
	CONVERTER(X, 8)\
	CONVERTER(X, 16)
#define CONVERTER_VECTORS_SAT(X)\
	CONVERTER_SAT(X, 2)\
	CONVERTER_SAT(X, 3)\
	CONVERTER_SAT(X, 4)\
	CONVERTER_SAT(X, 8)\
	CONVERTER_SAT(X, 16)

CONVERTER_VECTORS(char)
CONVERTER_VECTORS(uchar)
CONVERTER_VECTORS(short)
CONVERTER_VECTORS(ushort)
CONVERTER_VECTORS(int)
CONVERTER_VECTORS(uint)
CONVERTER_VECTORS(long)
CONVERTER_VECTORS(ulong)
CONVERTER_VECTORS(float)
CONVERTER_VECTORS(double)

CONVERTER_VECTORS_SAT(char)
CONVERTER_VECTORS_SAT(uchar)
CONVERTER_VECTORS_SAT(short)
CONVERTER_VECTORS_SAT(ushort)
CONVERTER_VECTORS_SAT(int)
CONVERTER_VECTORS_SAT(uint)
CONVERTER_VECTORS_SAT(long)
CONVERTER_VECTORS_SAT(ulong)

#undef CONVERTER_VECTORS
#undef CONVERTER
#undef CONVERTER_VECTORS_SAT
#undef CONVERTER_SAT

#define AS_TYPE(Out)\
template<typename In>	static inline const __##Out &as_##Out(const In &t)	{	return reinterpret_cast<const __##Out&>(t);	}\
template<typename In>	static inline __##Out &as_##Out(In &t)	{	return reinterpret_cast<__##Out&>(t);	}

#define AS_BASE_TYPE(Out)\
AS_TYPE(Out)\
AS_TYPE(Out##2)\
AS_TYPE(Out##3)\
AS_TYPE(Out##4)\
AS_TYPE(Out##8)\
AS_TYPE(Out##16)

AS_BASE_TYPE(char)
AS_BASE_TYPE(short)
AS_BASE_TYPE(int)
AS_BASE_TYPE(long)
AS_BASE_TYPE(uchar)
AS_BASE_TYPE(ushort)
AS_BASE_TYPE(uint)
AS_BASE_TYPE(ulong)
AS_BASE_TYPE(float)
AS_BASE_TYPE(double)
AS_TYPE(size_t)

#undef AS_TYPE
#undef AS_BASE_TYPE

#endif

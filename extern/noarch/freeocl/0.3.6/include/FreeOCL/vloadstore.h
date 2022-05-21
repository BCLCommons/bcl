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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_VLOADSTORE_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_VLOADSTORE_H__
// Built-in vector load and store functions
// vloadn
template<typename S>
static inline const typename __vector_type<S, 2>::type &vload2(size_t offset, const S *p)
{	return *((typename __vector_type<S, 2>::type*)(p + offset * 2));	}

template<typename S>
static inline typename __vector_type<S, 3>::type vload3(size_t offset, const S *p)
{
	typename __vector_type<S, 3>::type ret;
	p += offset * 3;
	ret.v[0] = p[0];
	ret.v[1] = p[1];
	ret.v[2] = p[2];
	return ret;
}

template<typename S>
static inline const typename __vector_type<S, 4>::type &vload4(size_t offset, const S *p)
{	return *((typename __vector_type<S, 4>::type*)(p + offset * 4));	}

template<typename S>
static inline const typename __vector_type<S, 8>::type &vload8(size_t offset, const S *p)
{	return *((typename __vector_type<S, 8>::type*)(p + offset * 8));	}

template<typename S>
static inline const typename __vector_type<S, 16>::type &vload16(size_t offset, const S *p)
{	return *((typename __vector_type<S, 16>::type*)(p + offset * 16));	}

// vstoren
template<typename S>
static inline void vstore2(const typename __vector_type<S, 2>::type &v, size_t offset, S *p)
{	*(typename __vector_type<S, 2>::type*)(p + offset * 2) = v;	}

template<typename S>
static inline void vstore3(const typename __vector_type<S, 3>::type &v, size_t offset, S *p)
{
	p += offset * 3;
	p[0] = v.v[0];
	p[1] = v.v[1];
	p[2] = v.v[2];
}

template<typename S>
static inline void vstore4(const typename __vector_type<S, 4>::type &v, size_t offset, S *p)
{	*(typename __vector_type<S, 4>::type*)(p + offset * 4) = v;	}

template<typename S>
static inline void vstore8(const typename __vector_type<S, 8>::type &v, size_t offset, S *p)
{	*(typename __vector_type<S, 8>::type*)(p + offset * 8) = v;	}

template<typename S>
static inline void vstore16(const typename __vector_type<S, 16>::type &v, size_t offset, S *p)
{	*(typename __vector_type<S, 16>::type*)(p + offset * 16) = v;	}

// vload_half
static inline float vload_half(size_t offset, const half *p)
{
	return p[offset];
}

// vstore_half_*
static inline void vstore_half(float data, size_t offset, half *p)
{
	p[offset] = half::from_float(data);
}

static inline void vstore_half_rte(float data, size_t offset, half *p)
{
	p[offset] = half::from_float(data);
}

#endif

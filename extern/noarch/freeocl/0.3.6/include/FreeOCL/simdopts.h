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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_SIMDOPTS_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_SIMDOPTS_H__

#ifdef __SSE__
#include <x86intrin.h>

static inline __float4 operator+(const __float4 &lhs, const __float4 &rhs)
{
	const __m128 s = _mm_add_ps(reinterpret_cast<const __m128&>(lhs), reinterpret_cast<const __m128&>(rhs));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 operator-(const __float4 &lhs, const __float4 &rhs)
{
	const __m128 s = _mm_sub_ps(reinterpret_cast<const __m128&>(lhs), reinterpret_cast<const __m128&>(rhs));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 operator*(const __float4 &lhs, const __float4 &rhs)
{
	const __m128 s = _mm_mul_ps(reinterpret_cast<const __m128&>(lhs), reinterpret_cast<const __m128&>(rhs));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 operator/(const __float4 &lhs, const __float4 &rhs)
{
	const __m128 s = _mm_div_ps(reinterpret_cast<const __m128&>(lhs), reinterpret_cast<const __m128&>(rhs));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 native_sqrt(const __float4 &v)
{
	const __m128 s = _mm_sqrt_ps(reinterpret_cast<const __m128&>(v));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 native_recip(const __float4 &v)
{
	const __m128 s = _mm_rcp_ps(reinterpret_cast<const __m128&>(v));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 native_rsqrt(const __float4 &v)
{
	const __m128 s = _mm_rsqrt_ps(reinterpret_cast<const __m128&>(v));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 native_divide(const __float4 &lhs, const __float4 &rhs)
{
	const __m128 s = _mm_div_ps(reinterpret_cast<const __m128&>(lhs), reinterpret_cast<const __m128&>(rhs));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 clamp(const __float4 &v, __float l, __float u)
{
	const __m128 s = _mm_min_ps(_mm_set1_ps(u), _mm_max_ps(reinterpret_cast<const __m128&>(v), _mm_set1_ps(l)));
	return reinterpret_cast<const __float4&>(s);
}

static inline __float4 clamp(const __float4 &v, const __float4 &l, const __float4 &u)
{
	const __m128 s = _mm_min_ps(reinterpret_cast<const __m128&>(u), _mm_max_ps(reinterpret_cast<const __m128&>(v), reinterpret_cast<const __m128&>(l)));
	return reinterpret_cast<const __float4&>(s);
}

#endif

#endif

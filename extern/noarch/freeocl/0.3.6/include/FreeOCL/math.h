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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_MATH_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_MATH_H__
#include "vectors.h"

#include <math.h>

// Built-in math functions
// for double
static inline __double acospi(__double x)	{	return acos(x) / M_PI;	}
static inline __double asinpi(__double x)	{	return asin(x) / M_PI;	}
static inline __double atanpi(__double x)	{	return atan(x) / M_PI;	}
static inline __double atan2pi(__double y, __double x)	{	return atan2(y, x) / M_PI;	}
static inline __double cospi(__double x)	{	return cos(M_PI * x);	}
static inline __double fract(__double x, __double *iptr)
{
	const __double _x = floor(x);
	*iptr = _x;
	return fmin(x - _x, 0x1.fffffep-1f);
}
static inline __double mad(__double a, __double b, __double c)	{	return a * b + c;	}
static inline __double maxmag(__double x, __double y)	{	return fmax(fabs(x), fabs(y));	}
static inline __double minmag(__double x, __double y)	{	return fmin(fabs(x), fabs(y));	}
static inline __double pown(__double x, __int n)
{
	__double _m = 1.0;
	if (n < 0)
	{
		x = 1.0 / x;
		n = -n;
	}
	for( ; n > 0 ; n >>= 1)
	{
		if (n & 1)
			_m *= x;
	}
	return _m;
}
static inline __double powr(__double x, __double y)	{	return pow(x, y);	}
static inline __double rootn(__double x, __int y)	{	return pow(x, 1.0 / y);	}
static inline __double rsqrt(__double x)	{	return 1.0 / sqrt(x);	}
static inline __double sincos(__double x, __double *cosval)
{
#ifdef _GNU_SOURCE
	__double _sin;
	sincos(x, &_sin, cosval);
	return _sin;
#else
	*cosval = cos(x);
	return sin(x);
#endif
}
static inline __double sinpi(__double x)	{	return sin(M_PI * x);	}
static inline __double tanpi(__double x)	{	return tan(M_PI * x);	}

// for float
static inline __float acos(__float x)	{	return acosf(x);	}
static inline __float acosh(__float x)	{	return acoshf(x);	}
static inline __float acospi(__float x)	{	return acos(x) / M_PI_F;	}
static inline __float asin(__float x)	{	return asinf(x);	}
static inline __float asinh(__float x)	{	return asinhf(x);	}
static inline __float asinpi(__float x)	{	return asin(x) / M_PI_F;	}
static inline __float atan(__float x)	{	return atanf(x);	}
static inline __float atan2(__float y, __float x)	{	return atan2f(y, x);	}
static inline __float atanh(__float x)	{	return atanhf(x);	}
static inline __float atanpi(__float x)	{	return atan(x) / M_PI_F;	}
static inline __float atan2pi(__float y, __float x)	{	return atan2(y, x) / M_PI_F;	}
static inline __float cbrt(__float x)	{	return cbrtf(x);	}
static inline __float ceil(__float x)	{	return ceilf(x);	}
static inline __float copysign(__float x, __float y)	{	return copysignf(x, y);	}
static inline __float cos(__float x)	{	return cosf(x);	}
static inline __float cosh(__float x)	{	return coshf(x);	}
static inline __float cospi(__float x)	{	return cos(M_PI_F * x);	}
static inline __float erfc(__float x)	{	return erfcf(x);	}
static inline __float erf(__float x)	{	return erff(x);	}
static inline __float exp(__float x)	{	return expf(x);	}
static inline __float exp2(__float x)	{	return exp2f(x);	}
static inline __float exp10(__float x)
{
#ifdef _GNU_SOURCE
	return exp10f(x);
#else
	return exp(x * 2.3025851f);
#endif
}
static inline __float expm1(__float x)	{	return expm1f(x);	}
static inline __float fabs(__float x)	{	return fabsf(x);	}
static inline __float fdim(__float x, __float y)	{	return fdimf(x, y);	}
static inline __float floor(__float x)	{	return floorf(x);	}
static inline __float fma(__float a, __float b, __float c)	{	return fmaf(a, b, c);	}
static inline __float fmax(__float x, __float y)	{	return fmaxf(x, y);	}
static inline __float fmin(__float x, __float y)	{	return fminf(x, y);	}
static inline __float fmod(__float x, __float y)	{	return fmodf(x, y);	}
static inline __float fract(__float x, __float *iptr)
{
	const __float _x = floor(x);
	*iptr = _x;
	return fmin(x - _x, 0x1.fffffep-1f);
}
static inline __float frexp(__float x, __int *exp)	{	return frexpf(x, exp);	}
static inline __float hypot(__float x, __float y)	{	return hypotf(x, y);	}
static inline __float ldexp(__float x, __int n)	{	return ldexpf(x, n);	}
static inline __float lgamma(__float x)	{	return lgammaf(x);	}
static inline __float lgamma_r(__float x, __int *signp)	{	return lgammaf_r(x, signp);	}
static inline __float log(__float x)	{	return logf(x);	}
static inline __float log2(__float x)	{	return log2f(x);	}
static inline __float log10(__float x)	{	return log10f(x);	}
static inline __float log1p(__float x)	{	return log1pf(x);	}
static inline __float logb(__float x)	{	return logbf(x);	}
static inline __float mad(__float a, __float b, __float c)	{	return a * b + c;	}
static inline __float maxmag(__float x, __float y)	{	return fmax(fabs(x), fabs(y));	}
static inline __float minmag(__float x, __float y)	{	return fmin(fabs(x), fabs(y));	}
static inline __float modf(__float x, __float *iptr)	{	return modff(x, iptr);	}
static inline __float nan(__uint nancode)	{	return reinterpret_cast<const __float&>(nancode |= 0x7F800000U);	}
static inline __float nextafter(__float x, __float y)	{	return nextafterf(x, y);	}
static inline __float pow(__float x, __float y)	{	return powf(x, y);	}
static inline __float pown(__float x, __int n)
{
	__float _m = 1.0f;
	if (n < 0)
	{
		x = 1.0f / x;
		n = -n;
	}
	for( ; n > 0 ; n >>= 1)
	{
		if (n & 1)
			_m *= x;
	}
	return _m;
}
static inline __float powr(__float x, __float y)	{	return powf(x, y);	}
static inline __float remainder(__float x, __float y)	{	return remainderf(x, y);	}
static inline __float remquo(__float x, __float y, __int *quo)	{	return remquof(x, y, quo);	}
static inline __float rint(__float x)	{	return rintf(x);	}
static inline __float rootn(__float x, __int y)	{	return pow(x, 1.0f / y);	}
static inline __float round(__float x)	{	return roundf(x);	}
static inline __float rsqrt(__float x)	{	return 1.0f / sqrtf(x);	}
static inline __float sin(__float x)	{	return sinf(x);	}
static inline __float sincos(__float x, __float *cosval)
{
#ifdef _GNU_SOURCE
	__float _sin;
	sincosf(x, &_sin, cosval);
	return _sin;
#else
	*cosval = cos(x);
	return sin(x);
#endif
}
static inline __float sinh(__float x)	{	return sinhf(x);	}
static inline __float sinpi(__float x)	{	return sin(M_PI_F * x);	}
static inline __float sqrt(__float x)	{	return sqrtf(x);	}
static inline __float tan(__float x)	{	return tanf(x);	}
static inline __float tanh(__float x)	{	return tanhf(x);	}
static inline __float tanpi(__float x)	{	return tan(M_PI_F * x);	}
static inline __float tgamma(__float x)	{	return tgammaf(x);	}
static inline __float trunc(__float x)	{	return truncf(x);	}

// half_ versions
static inline __float half_cos(__float x)	{	return cos(x);	}
static inline __float half_divide(__float x, __float y)	{	return x / y;	}
static inline __float half_exp(__float x)	{	return exp(x);	}
static inline __float half_exp2(__float x)	{	return exp2(x);	}
static inline __float half_exp10(__float x)	{	return exp10(x);	}
static inline __float half_log(__float x)	{	return log(x);	}
static inline __float half_log2(__float x)	{	return log2(x);	}
static inline __float half_log10(__float x)	{	return log10(x);	}
static inline __float half_powr(__float x, __float y)	{	return powr(x, y);	}
static inline __float half_recip(__float x)	{	return 1.0f / x;	}
static inline __float half_rsqrt(__float x)	{	return 1.0f / sqrt(x);	}
static inline __float half_sin(__float x)	{	return sin(x);	}
static inline __float half_sqrt(__float x)	{	return sqrt(x);	}
static inline __float half_tan(__float x)	{	return tan(x);	}

// native_ versions
//static inline __float native_cos(__float x)	{	return cos(x);	}
static inline __float native_divide(__float x, __float y)	{	return x / y;	}
//static inline __float native_exp(__float x)	{	return exp(x);	}
static inline __float native_exp2(__float x)	{	return exp2(x);	}
static inline __float native_exp10(__float x)	{	return exp10(x);	}
//static inline __float native_log(__float x)	{	return log(x);	}
static inline __float native_log2(__float x)	{	return log2(x);	}
static inline __float native_log10(__float x)	{	return log10(x);	}
//static inline __float native_powr(__float x, __float y)	{	return powr(x, y);	}
static inline __float native_recip(__float x)	{	return 1.0f / x;	}
static inline __float native_rsqrt(__float x)	{	return 1.0f / sqrt(x);	}
//static inline __float native_sin(__float x)	{	return sin(x);	}
static inline __float native_sqrt(__float x)	{	return sqrt(x);	}
//static inline __float native_tan(__float x)	{	return tan(x);	}

static inline __float native_sin(const __float x)
{
	__float sign, y;
	if (x < 0.0f)
	{
		sign = -1.0f;
		y = -x;
	}
	else
	{
		sign = 1.0f;
		y = x;
	}
	__int q = y * __float(1.0 / M_PI_2);
	y -= q * __float(M_PI_2);
	if (q & 2)
		sign = -sign;
	if (q & 1)
		y = __float(M_PI_2) - y;

	const __float y2 = y * y;
	__float z;
	z = __float(1.0 / 362880.0);
	z = z * y2 - __float(1.0 / 5040.0);
	z = z * y2 + __float(1.0 / 120.0);
	z = z * y2 - __float(1.0 / 6.0);
	z = z * y2 + 1.0f;
	z *= y;
	return sign * z;
}

static inline __float native_cos(const __float x)
{
	__float sign = 1.0f, y;
	if (x < 0.0f)
		y = -x;
	else
		y = x;
	__int q = y * __float(1.0 / M_PI_2);
	y -= q * __float(M_PI_2);
	q &= 3;
	if (q == 1 || q == 2)
		sign = -sign;
	if (q & 1)
		y = __float(M_PI_2) - y;

	y *= y;
	__float z;
	z = __float(1.0 / 40320.0);
	z = z * y - __float(1.0 / 720.0);
	z = z * y + __float(1.0 / 24.0);
	z = z * y - __float(1.0 / 2.0);
	z = z * y + 1.0f;
	return sign * z;
}

static inline __float native_ldexp(const __float x, const __int exp)
{
	const __uint ix = reinterpret_cast<const __uint &>(x);
	const __uint e = (ix & 0x7F800000U) + (exp << 23);
	const __uint f = (ix & ~0x7F800000U) | e;
	return reinterpret_cast<const __float&>(f);
}

static inline __float native_exp(const __float x)
{
//	if (x > 88.7228f)
//		return __float(INFINITY);
//	if (x < -87.3365f)
//		return 0.0f;
	__float y = x * __float(1.0 / M_LN2);
	const __int j = y;
	y = fmodf(y, 1.0f);
	y *= __float(M_LN2);
	__float z;
	z = __float(1.0 / 720.0);
	z = z * y  + __float(1.0 / 120.0);
	z = z * y  + __float(1.0 / 24.0);
	z = z * y  + __float(1.0 / 6.0);
	z = z * y  + __float(1.0 / 2.0);
	z = z * y  + 1.0f;
	z = z * y  + 1.0f;
	return native_ldexp(z, j);
}

static inline __float native_log(const __float x)
{
	__int e;
	const __float y = frexpf(x, &e) - 1.0f;
	__float z;
	z = -__float(1.0 / 10.0);
	z = z * y  + __float(1.0 / 9.0);
	z = z * y  - __float(1.0 / 8.0);
	z = z * y  + __float(1.0 / 7.0);
	z = z * y  - __float(1.0 / 6.0);
	z = z * y  + __float(1.0 / 5.0);
	z = z * y  - __float(1.0 / 4.0);
	z = z * y  + __float(1.0 / 3.0);
	z = z * y  - __float(1.0 / 2.0);
	z = z * y  + 1.0f;
	z = z * y;
	return z + e * __float(M_LN2);
}

static inline __float native_tan(const __float x)
{
	return native_sin(x) / native_cos(x);
}

static inline __float native_powr(const __float x, const __float y)
{
	return native_exp(y * native_log(x));
}

// Vector versions
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(acos)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(acosh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(acospi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(asin)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(asinh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(asinpi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(atan)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(atan2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(atanh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(atanpi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(atan2pi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(cbrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(ceil)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(copysign)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(cos)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(cosh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(cospi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(erfc)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(erf)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(exp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(exp2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(exp10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(expm1)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(fabs)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(fdim)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(floor)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3(fma)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(fmax)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(fmin)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(fmod)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PB(fract)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PI(frexp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(hypot)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2I(ldexp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(lgamma)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PI(lgamma_r)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(log)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(log2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(log10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(log1p)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(logb)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3(mad)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(maxmag)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(minmag)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PB(modf)
static inline __float2 nan(__uint2 nancode)	{	return __float2::make(nan(nancode.v[0]), nan(nancode.v[1]));	}
static inline __float3 nan(__uint3 nancode)	{	return __float3::make(nan(nancode.v[0]), nan(nancode.v[1]), nan(nancode.v[2]));	}
static inline __float4 nan(__uint4 nancode)	{	return __float4::make(nan(nancode.v[0]), nan(nancode.v[1]), nan(nancode.v[2]), nan(nancode.v[3]));	}
static inline __float8 nan(__uint8 nancode)	{	return __float8::make(nan(nancode.v[0]), nan(nancode.v[1]), nan(nancode.v[2]), nan(nancode.v[3]),
																  nan(nancode.v[4]), nan(nancode.v[5]), nan(nancode.v[6]), nan(nancode.v[7]));	}
static inline __float16 nan(__uint16 nancode)	{	return __float16::make(nan(nancode.v[0]), nan(nancode.v[1]), nan(nancode.v[2]), nan(nancode.v[3]),
																   nan(nancode.v[4]), nan(nancode.v[5]), nan(nancode.v[6]), nan(nancode.v[7]),
																   nan(nancode.v[8]), nan(nancode.v[9]), nan(nancode.v[10]), nan(nancode.v[11]),
																   nan(nancode.v[12]), nan(nancode.v[13]), nan(nancode.v[14]), nan(nancode.v[15]));	}
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(nextafter)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(pow)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2I(pown)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(powr)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(remainder)
template<class V> static inline CHECK(V) remquo(const V &x, const V &y, const typename __vector_type<__int, __vector<V>::components>::type &z)
{
	V ret;
	for(size_t i = 0 ; i < __vector<V>::components ; ++i)
		ret.v[i] = F(x.v[i], y.v[i], (__int*)z + i);
	return ret;
}
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(rint)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2I(rootn)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(round)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(rsqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(sin)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2PB(sincos)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(sinh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(sinpi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(sqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(tan)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(tanh)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(tanpi)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(tgamma)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(trunc)

// half_ versions
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_cos)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(half_divide)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_exp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_exp2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_exp10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_log)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_log2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_log10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(half_powr)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_recip)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_rsqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_sin)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_sqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(half_tan)

// native_ versions
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_cos)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(native_divide)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_exp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_exp2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_exp10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_log)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_log2)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_log10)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(native_powr)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_recip)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_rsqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_sin)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_sqrt)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(native_tan)

#endif

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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_COMMON_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_COMMON_H__
#include "vectors.h"

// Built-in common functions
// for __float
static inline __float clamp(__float x, __float minval, __float maxval)	{	return fmin(fmax(x, minval), maxval);	}
static inline __float degrees(__float radians)	{	return radians * 57.295779513082f;	}
static inline __float max(__float x, __float y)	{	return x < y ? y : x;	}
static inline __float min(__float x, __float y)	{	return y < x ? y : x;	}
static inline __float mix(__float x, __float y, __float a)	{	return x + (y - x) * a;	}
static inline __float radians(__float degrees)	{	return degrees * 0.017453292519943f;	}
static inline __float step(__float edge, __float x)	{	return x < edge ? 0.0f : 1.0f;	}
static inline __float smoothstep(__float edge0, __float edge1, __float x)
{
	__float t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
	return t * t * (3.0f - 2.0f * t);
}
static inline __float sign(__float x)
{
	return x > 0.0f ? 1.0f : x < 0.0f ? -1.0f : isnan(x) ? 0.0f : x;
}

// for __double
static inline __double clamp(__double x, __double minval, __double maxval)	{	return fmin(fmax(x, minval), maxval);	}
static inline __double degrees(__double radians)	{	return radians * 57.295779513082;	}
static inline __double max(__double x, __double y)	{	return x < y ? y : x;	}
static inline __double min(__double x, __double y)	{	return y < x ? y : x;	}
static inline __double mix(__double x, __double y, __double a)	{	return x + (y - x) * a;	}
static inline __double radians(__double degrees)	{	return degrees * 0.017453292519943;	}
static inline __double step(__double edge, __double x)	{	return x < edge ? 0.0 : 1.0;	}
static inline __double smoothstep(__double edge0, __double edge1, __double x)
{
	__double t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
	return t * t * (3.0 - 2.0 * t);
}
static inline __double sign(__double x)
{
	return x > 0.0 ? 1.0 : x < 0.0 ? -1.0 : isnan(x) ? 0.0 : x;
}

VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3(clamp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3SS(clamp)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(degrees)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(max)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2S(max)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(min)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2S(min)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3(mix)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3_S(mix)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(radians)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2(step)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION2S_(step)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3(smoothstep)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION3SS_(smoothstep)
VECTOR_IMPLEMENTATION_FROM_SCALAR_IMPLEMENTATION1(sign)

#endif

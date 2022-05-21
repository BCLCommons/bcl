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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_IMGREADWRITE_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_IMGREADWRITE_H__

typedef __uint sampler_t;

#define CLK_NORMALIZED_COORDS_FALSE		0x00000000U
#define CLK_NORMALIZED_COORDS_TRUE		0x01000000U

#define CLK_ADDRESS_NONE				0x00000000U
#define CLK_ADDRESS_CLAMP_TO_EDGE		0x00000001U
#define CLK_ADDRESS_REPEAT				0x00000002U
#define CLK_ADDRESS_CLAMP				0x00000003U
#define CLK_ADDRESS_MIRRORED_REPEAT		0x00000004U

#define CLK_FILTER_NEAREST				0x00000000U
#define CLK_FILTER_LINEAR				0x00010000U

/* cl_channel_order */
#define CLK_R			0x10B0
#define CLK_A			0x10B1
#define CLK_RG			0x10B2
#define CLK_RA			0x10B3
#define CLK_RGB			0x10B4
#define CLK_RGBA		0x10B5
#define CLK_BGRA		0x10B6
#define CLK_ARGB		0x10B7
#define CLK_INTENSITY	0x10B8
#define CLK_LUMINANCE	0x10B9
#define CLK_Rx			0x10BA
#define CLK_RGx			0x10BB
#define CLK_RGBx		0x10BC

/* cl_channel_type */
#define CLK_SNORM_INT8			0x10D0
#define CLK_SNORM_INT16			0x10D1
#define CLK_UNORM_INT8			0x10D2
#define CLK_UNORM_INT16			0x10D3
#define CLK_UNORM_SHORT_565		0x10D4
#define CLK_UNORM_SHORT_555		0x10D5
#define CLK_UNORM_INT_101010	0x10D6
#define CLK_SIGNED_INT8			0x10D7
#define CLK_SIGNED_INT16		0x10D8
#define CLK_SIGNED_INT32		0x10D9
#define CLK_UNSIGNED_INT8		0x10DA
#define CLK_UNSIGNED_INT16		0x10DB
#define CLK_UNSIGNED_INT32		0x10DC
#define CLK_HALF_FLOAT			0x10DD
#define CLK_FLOAT				0x10DE

static inline __float __third(__float v)	{	return v * 0.3333333333f;	}
static inline __int __third(__int v)	{	return v / 3;	}
static inline __uint __third(__uint v)	{	return v / 3U;	}
static inline __float __quarter(__float v)	{	return v * 0.25f;	}
static inline __int __quarter(__int v)	{	return v >> 2;	}
static inline __uint __quarter(__uint v)	{	return v >> 2;	}

template<class I, class T>
static inline T __map_channels_for_reading(const I &image, const T &v)
{
	T color;
	switch(image.channel_order)
	{
	case CLK_R:
	case CLK_Rx:
		color = T::make(v.v[0],0,0,0);
		break;
	case CLK_A:
		color = T::make(0,0,0,v.v[0]);
		break;
	case CLK_INTENSITY:
		color = T::make(v.v[0]);
		break;
	case CLK_LUMINANCE:
		color = T::make(v.v[0],v.v[0],v.v[0],0);
		break;

	case CLK_RG:
	case CLK_RGx:
		color = T::make(v.v[0],v.v[1],0,0);
		break;
	case CLK_RA:
		color = T::make(v.v[0],0,0,v.v[1]);
		break;

	case CLK_RGB:
	case CLK_RGBx:
		color = T::make(v.v[0],v.v[1],v.v[2],0);
		break;

	case CLK_RGBA:
		color = v;
		break;
	case CLK_BGRA:
		color = T::make(v.v[2],v.v[1],v.v[0],v.v[3]);
		break;
	case CLK_ARGB:
		color = T::make(v.v[1],v.v[2],v.v[3],v.v[0]);
		break;
	}
	return color;
}

template<typename I, class T>
static inline T __map_channels_for_writing(const I &image, const T &color, __size_t &nb_chan)
{
	T v;
	switch(image.channel_order)
	{
	case CLK_R:
		v.v[0] = color.v[0];
		nb_chan = 1;
		break;
	case CLK_A:
		v.v[0] = color.v[3];
		nb_chan = 1;
		break;
	case CLK_INTENSITY:
		v.v[0] = __quarter(color.v[0] + color.v[1] + color.v[2] + color.v[3]);
		nb_chan = 1;
		break;
	case CLK_LUMINANCE:
		v.v[0] = __third(color.v[0] + color.v[1] + color.v[2]);
		nb_chan = 1;
		break;
	case CLK_Rx:
		v.v[0] = color.v[0];
		nb_chan = 1;
		break;

	case CLK_RG:
		v.v[0] = color.v[0];
		v.v[1] = color.v[1];
		nb_chan = 2;
		break;
	case CLK_RA:
		v.v[0] = color.v[0];
		v.v[1] = v.v[3];
		nb_chan = 2;
		break;
	case CLK_RGx:
		v.v[0] = color.v[0];
		v.v[1] = color.v[1];
		nb_chan = 2;
		break;

	case CLK_RGB:
		v.v[0] = color.v[0];
		v.v[1] = color.v[1];
		v.v[2] = color.v[2];
		nb_chan = 3;
		break;
	case CLK_RGBx:
		v.v[0] = color.v[0];
		v.v[1] = color.v[1];
		v.v[2] = color.v[2];
		nb_chan = 3;
		break;

	case CLK_RGBA:
		v = color;
		nb_chan = 4;
		break;
	case CLK_BGRA:
		v.v[0] = color.v[2];
		v.v[1] = color.v[1];
		v.v[2] = color.v[0];
		v.v[3] = color.v[3];
		nb_chan = 4;
		break;
	case CLK_ARGB:
		v.v[0] = color.v[3];
		v.v[1] = color.v[0];
		v.v[2] = color.v[1];
		v.v[3] = color.v[2];
		nb_chan = 4;
		break;
	}
	return v;
}

#include "image1d_t.h"
#include "image1d_buffer_t.h"
#include "image1d_array_t.h"
#include "image2d_t.h"
#include "image2d_array_t.h"
#include "image3d_t.h"

#endif

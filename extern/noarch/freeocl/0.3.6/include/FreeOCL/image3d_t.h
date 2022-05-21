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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_IMAGE3D_T_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_IMAGE3D_T_H__

struct image3d_t
{
	__uint channel_order;
	__uint channel_data_type;
	__size_t width, height, depth;
	__size_t row_pitch, slice_pitch;
	__size_t element_size;
	void *data;
};

static inline __int4 __address_mode(const image3d_t &image, const sampler_t sampler, const __int4 &coord)
{
	__int4 __coord;
	switch(sampler & 0xFF)
	{
	case CLK_ADDRESS_NONE:
		__coord = coord;
		break;
	case CLK_ADDRESS_CLAMP_TO_EDGE:
		__coord = clamp(coord, __int4::make(0), __int4::make(image.width - 1, image.height - 1, image.depth - 1, 0));
		break;
	case CLK_ADDRESS_CLAMP:
		__coord = clamp(coord, __int4::make(-1), __int4::make(image.width, image.height, image.depth, 0));
		break;
	case CLK_ADDRESS_REPEAT:
		__coord = coord % __int4::make(image.width, image.height, image.depth, 1);
		break;
	case CLK_ADDRESS_MIRRORED_REPEAT:
		__coord = coord % __int4::make(image.width << 1, image.height << 1, image.depth << 1, 1);
		__coord = abs(__coord - __int4::make(image.width, image.height, image.depth, 0));
		break;
	}
	return __coord;
}

static inline __float4 read_imagef(const image3d_t &image, sampler_t sampler, const __int4 &coord)
{
	const __int4 __coord = __address_mode(image, sampler, coord);

	const __char * const ptr = (const __char*)image.data
							 + __coord.get<0>() * image.element_size
							 + __coord.get<1>() * image.row_pitch
							 + __coord.get<2>() * image.slice_pitch;
	__size_t nb_chan;
	switch(image.channel_order)
	{
	case CLK_R:			nb_chan = 1;	break;
	case CLK_A:			nb_chan = 1;	break;
	case CLK_INTENSITY:	nb_chan = 1;	break;
	case CLK_LUMINANCE:	nb_chan = 1;	break;
	case CLK_Rx:		nb_chan = 1;	break;

	case CLK_RG:		nb_chan = 2;	break;
	case CLK_RA:		nb_chan = 2;	break;
	case CLK_RGx:		nb_chan = 2;	break;

	case CLK_RGB:		nb_chan = 3;	break;
	case CLK_RGBx:		nb_chan = 3;	break;

	case CLK_RGBA:		nb_chan = 4;	break;
	case CLK_BGRA:		nb_chan = 4;	break;
	case CLK_ARGB:		nb_chan = 4;	break;
	}

	__float4 v;
	switch(image.channel_data_type)
	{
	case CLK_UNORM_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __uchar*)ptr)[i] * (1.0f / 255.0f);
		break;
	case CLK_UNORM_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __ushort*)ptr)[i] * (1.0f / 65535.0f);
		break;
	case CLK_UNORM_SHORT_555:
		{
			const __uint s = *(const __ushort*)ptr;
			v.get<0>() = ((s >> 10) & 31U) * (1.0f / 31.0f);
			v.get<1>() = ((s >> 5) & 31U) * (1.0f / 31.0f);
			v.get<2>() = (s & 31U) * (1.0f / 31.0f);
		}
		break;
	case CLK_UNORM_SHORT_565:
		{
			const __uint s = *(const __ushort*)ptr;
			v.get<0>() = ((s >> 11) & 31U) * (1.0f / 31.0f);
			v.get<1>() = ((s >> 5) & 63U) * (1.0f / 63.0f);
			v.get<2>() = (s & 31U) * (1.0f / 31.0f);
		}
		break;
	case CLK_UNORM_INT_101010:
		{
			const __uint s = *(const __uint*)ptr;
			v.get<0>() = ((s >> 20) & 1023U) * (1.0f / 1023.0f);
			v.get<1>() = ((s >> 10) & 1023U) * (1.0f / 1023.0f);
			v.get<2>() = (s & 1023U) * (1.0f / 1023.0f);
		}
		break;
	case CLK_SNORM_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __char*)ptr)[i] * (1.0f / 127.0f);
		break;
	case CLK_SNORM_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __short*)ptr)[i] * (1.0f / 32767.0f);
		break;
	case CLK_HALF_FLOAT:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __half*)ptr)[i];
		break;
	case CLK_FLOAT:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __float*)ptr)[i];
		break;
	}

	return __map_channels_for_reading(image, v);
}
static inline __float4 read_imagef(const image3d_t &image, sampler_t sampler, const __float4 &coord)
{
	__float4 rcoord = coord;
	if (sampler & CLK_NORMALIZED_COORDS_TRUE)
		rcoord *= __float4::make(image.width, image.height, image.depth, 0.0f);
	if (!(sampler & CLK_FILTER_LINEAR))
		return read_imagef(image, sampler, convert_int4(rcoord));
	const __int4 c0 = convert_int4(rcoord - __float4::make(0.5f, 0.5f, 0.5f, 0.f));
	const __float4 dp = rcoord - convert_float4(c0);
	__float4 v0 = read_imagef(image, sampler, c0);
	__float4 v1 = read_imagef(image, sampler, c0 + __int4::make(1,0,0,0));
	v0 += dp.get<0>() * (v1 - v0);
	v1 = read_imagef(image, sampler, c0 + __int4::make(0,1,0,0));
	__float4 v2 = read_imagef(image, sampler, c0 + __int4::make(1,1,0,0));
	v1 += dp.get<0>() * (v2 - v1);
	v0 += dp.get<1>() * (v1 - v0);

	v1 = read_imagef(image, sampler, c0 + __int4::make(0,0,1,0));
	v2 = read_imagef(image, sampler, c0 + __int4::make(1,0,1,0));
	v1 += dp.get<0>() * (v2 - v1);
	v2 = read_imagef(image, sampler, c0 + __int4::make(0,1,1,0));
	const __float4 v3 = read_imagef(image, sampler, c0 + __int4::make(1,1,1,0));
	v2 += dp.get<0>() * (v3 - v2);
	v1 += dp.get<1>() * (v2 - v1);

	return v0 + dp.get<2>() * (v1 - v0);
}
static inline __float4 read_imagef(const image3d_t &image, const __int4 &coord)
{
	return read_imagef(image, CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE, coord);
}

static inline __int4 read_imagei(const image3d_t &image, sampler_t sampler, const __int4 &coord)
{
	const __int4 __coord = __address_mode(image, sampler, coord);

	const __char * const ptr = (const __char*)image.data
							 + __coord.get<0>() * image.element_size
							 + __coord.get<1>() * image.row_pitch
							 + __coord.get<2>() * image.slice_pitch;
	__size_t nb_chan;
	switch(image.channel_order)
	{
	case CLK_R:			nb_chan = 1;	break;
	case CLK_A:			nb_chan = 1;	break;
	case CLK_INTENSITY:	nb_chan = 1;	break;
	case CLK_LUMINANCE:	nb_chan = 1;	break;
	case CLK_Rx:		nb_chan = 1;	break;

	case CLK_RG:		nb_chan = 2;	break;
	case CLK_RA:		nb_chan = 2;	break;
	case CLK_RGx:		nb_chan = 2;	break;

	case CLK_RGB:		nb_chan = 3;	break;
	case CLK_RGBx:		nb_chan = 3;	break;

	case CLK_RGBA:		nb_chan = 4;	break;
	case CLK_BGRA:		nb_chan = 4;	break;
	case CLK_ARGB:		nb_chan = 4;	break;
	}

	__int4 v;
	switch(image.channel_data_type)
	{
	case CLK_SIGNED_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __char*)ptr)[i];
		break;
	case CLK_SIGNED_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __short*)ptr)[i];
		break;
	case CLK_SIGNED_INT32:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __int*)ptr)[i];
		break;
	}

	return __map_channels_for_reading(image, v);
}
static inline __int4 read_imagei(const image3d_t &image, sampler_t sampler, const __float4 &coord)
{
	if (sampler & CLK_NORMALIZED_COORDS_TRUE)
		return read_imagei(image, sampler, convert_int4(coord * __float4::make(image.width, image.height, image.depth, 0)));
	return read_imagei(image, sampler, convert_int4(coord));
}
static inline __int4 read_imagei(const image3d_t &image, const __int4 &coord)
{
	return read_imagei(image, CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE, coord);
}

static inline __uint4 read_imageui(const image3d_t &image, sampler_t sampler, const __int4 &coord)
{
	const __int4 __coord = __address_mode(image, sampler, coord);

	const __char * const ptr = (const __char*)image.data
							 + __coord.get<0>() * image.element_size
							 + __coord.get<1>() * image.row_pitch
							 + __coord.get<2>() * image.slice_pitch;
	__size_t nb_chan;
	switch(image.channel_order)
	{
	case CLK_R:			nb_chan = 1;	break;
	case CLK_A:			nb_chan = 1;	break;
	case CLK_INTENSITY:	nb_chan = 1;	break;
	case CLK_LUMINANCE:	nb_chan = 1;	break;
	case CLK_Rx:		nb_chan = 1;	break;

	case CLK_RG:		nb_chan = 2;	break;
	case CLK_RA:		nb_chan = 2;	break;
	case CLK_RGx:		nb_chan = 2;	break;

	case CLK_RGB:		nb_chan = 3;	break;
	case CLK_RGBx:		nb_chan = 3;	break;

	case CLK_RGBA:		nb_chan = 4;	break;
	case CLK_BGRA:		nb_chan = 4;	break;
	case CLK_ARGB:		nb_chan = 4;	break;
	}

	__uint4 v;
	switch(image.channel_data_type)
	{
	case CLK_UNSIGNED_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __uchar*)ptr)[i];
		break;
	case CLK_UNSIGNED_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __ushort*)ptr)[i];
		break;
	case CLK_UNSIGNED_INT32:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			v.v[i] = ((const __uint*)ptr)[i];
		break;
	}

	return __map_channels_for_reading(image, v);
}
static inline __uint4 read_imageui(const image3d_t &image, sampler_t sampler, const __float4 &coord)
{
	if (sampler & CLK_NORMALIZED_COORDS_TRUE)
		return read_imageui(image, sampler, convert_int4(coord * __float4::make(image.width, image.height, image.depth, 0)));
	return read_imageui(image, sampler, convert_int4(coord));
}
static inline __uint4 read_imageui(const image3d_t &image, const __int4 &coord)
{
	return read_imageui(image, CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE, coord);
}

static inline void write_imagef(image3d_t &image, const __int4 &coord, const __float4 &color)
{
	__size_t nb_chan;
	const __float4 v = __map_channels_for_writing(image, color, nb_chan);
	__char *const ptr = (__char*)image.data
					  + image.row_pitch * coord.get<1>()
					  + coord.get<0>() * image.element_size
					  + coord.get<2>() * image.slice_pitch;
	switch(image.channel_data_type)
	{
	case CLK_SNORM_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__char*)ptr)[i] = clamp(v.v[i] * 127.0f, -127.0f, 127.0f);
		break;
	case CLK_SNORM_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__short*)ptr)[i] = clamp(v.v[i] * 32767.0f, -32767.0f, 32767.0f);
		break;
	case CLK_UNORM_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__uchar*)ptr)[i] = clamp(v.v[i] * 255.0f, 0.0f, 255.0f);
		break;
	case CLK_UNORM_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__ushort*)ptr)[i] = clamp(v.v[i] * 65535.0f, 0.0f, 65535.0f);
		break;
	case CLK_UNORM_SHORT_565:
		*((__ushort*)ptr) = (((__uint)clamp(v.get<0>() * 31.0f, 0.0f, 31.0f)) << 11)
						| (((__uint)clamp(v.get<1>() * 63.0f, 0.0f, 63.0f)) << 5)
						| ((__uint)clamp(v.get<2>() * 31.0f, 0.0f, 31.0f));
		break;
	case CLK_UNORM_SHORT_555:
		*((__ushort*)ptr) = (((__uint)clamp(v.get<0>() * 31.0f, 0.0f, 31.0f)) << 10)
						| (((__uint)clamp(v.get<1>() * 31.0f, 0.0f, 31.0f)) << 5)
						| ((__uint)clamp(v.get<2>() * 31.0f, 0.0f, 31.0f));
		break;
	case CLK_UNORM_INT_101010:
		*((__uint*)ptr) = (((__uint)clamp(v.get<0>() * 1023.0f, 0.0f, 1023.0f)) << 20)
						| (((__uint)clamp(v.get<1>() * 1023.0f, 0.0f, 1023.0f)) << 10)
						| ((__uint)clamp(v.get<2>() * 1023.0f, 0.0f, 1023.0f));
		break;
	case CLK_HALF_FLOAT:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__half*)ptr)[i] = __half::from_float(v.v[i]);
		break;
	case CLK_FLOAT:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__float*)ptr)[i] = v.v[i];
		break;
	}
}
static inline void write_imagei(image3d_t &image, const __int4 &coord, const __int4 &color)
{
	__size_t nb_chan;
	const __int4 v = __map_channels_for_writing(image, color, nb_chan);
	__char *const ptr = (__char*)image.data
					  + image.row_pitch * coord.get<1>()
					  + coord.get<0>() * image.element_size
					  + coord.get<2>() * image.slice_pitch;
	switch(image.channel_data_type)
	{
	case CLK_SIGNED_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__char*)ptr)[i] = clamp(v.v[i], -128, 127);
		break;
	case CLK_SIGNED_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__short*)ptr)[i] = clamp(v.v[i], -32768, 32767);
		break;
	case CLK_SIGNED_INT32:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__int*)ptr)[i] = v.v[i];
		break;
	}
}
static inline void write_imageui(image3d_t &image, const __int4 &coord, const __uint4 &color)
{
	__size_t nb_chan;
	const __uint4 v = __map_channels_for_writing(image, color, nb_chan);
	__char *const ptr = (__char*)image.data
					  + image.row_pitch * coord.get<1>()
					  + coord.get<0>() * image.element_size
					  + coord.get<2>() * image.slice_pitch;
	switch(image.channel_data_type)
	{
	case CLK_UNSIGNED_INT8:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__uchar*)ptr)[i] = clamp(v.v[i], 0, 255);
		break;
	case CLK_UNSIGNED_INT16:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__ushort*)ptr)[i] = clamp(v.v[i], 0, 65535);
		break;
	case CLK_SIGNED_INT32:
		for(__size_t i = 0 ; i < nb_chan ; ++i)
			((__uint*)ptr)[i] = v.v[i];
		break;
	}
}

static inline __int get_image_width(const image3d_t &image)	{	return image.width;	}
static inline __int get_image_height(const image3d_t &image)	{	return image.height;	}
static inline __int get_image_depth(const image3d_t &image)	{	return image.depth;	}
static inline __int get_image_channel_data_type(const image3d_t &image)	{	return image.channel_data_type;	}
static inline __int get_image_channel_order(const image3d_t &image)	{	return image.channel_order;	}

static inline __int4 get_image_dim(const image3d_t &image)
{
	return __int4::make(image.width, image.height, image.depth, 0);
}

#endif

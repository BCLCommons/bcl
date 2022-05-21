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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_WORKITEM_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_WORKITEM_H__
#include "FreeOCL/config.h"
#include <ucontext.h>

namespace FreeOCL
{
	static const void *args;
	static __thread char *local_memory;
	static __thread ucontext_t *scheduler;
	static __thread ucontext_t *threads;
	static __uint dim;
	static __size_t global_size[3];
	static __size_t global_offset[3];
	static __size_t local_size[3];
	static __thread __size_t group_id[3];
	static __size_t num_groups[3];
	static __thread __size_t thread_num;
}

// Built-in work-item functions
static inline __uint get_work_dim()
{
	return FreeOCL::dim;
}

static inline __size_t get_global_size(__uint dimindx)
{
	return dimindx < FreeOCL::dim ? FreeOCL::global_size[dimindx] : 1;
}

static inline __size_t get_local_id(__uint dimindx);
static inline __size_t get_global_id(__uint dimindx)
{
	return dimindx < FreeOCL::dim ? FreeOCL::group_id[dimindx] * FreeOCL::local_size[dimindx] + get_local_id(dimindx) + FreeOCL::global_offset[dimindx] : 0;
}

static inline __size_t get_local_size(__uint dimindx)
{
	return dimindx < FreeOCL::dim ? FreeOCL::local_size[dimindx] : 1;
}
static inline __size_t get_thread_num();
static inline __size_t get_local_id(__uint dimindx)
{
	switch(dimindx)
	{
	case 0:
		return get_thread_num() % FreeOCL::local_size[0];
	case 1:
		return (get_thread_num() / FreeOCL::local_size[0]) % FreeOCL::local_size[1];
	case 2:
		return (get_thread_num() / FreeOCL::local_size[0]) / FreeOCL::local_size[1];
	default:
		return -1;
	}
}

static inline __size_t get_num_groups(__uint dimindx)
{
	return dimindx < FreeOCL::dim ? FreeOCL::num_groups[dimindx] : 1;
}

static inline __size_t get_group_id(__uint dimindx)
{
	return FreeOCL::group_id[dimindx];
}

static inline __size_t get_global_offset(__uint dimindx)
{
	return dimindx < FreeOCL::dim ? FreeOCL::global_offset[dimindx] : 0;
}

static inline __size_t get_thread_num()
{
	return FreeOCL::thread_num;
}
#endif

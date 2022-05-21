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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_VMISC_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_VMISC_H__
// Built-in miscellaneous vector functions
// implement vec_step as a macro using templates and typeof
template<class T> struct vec_step { enum { value = 1 }; };
template<class T> struct vec_step<__vector_type<T, 2> > { enum { value = 2 }; };
template<class T> struct vec_step<__vector_type<T, 3> > { enum { value = 4 }; };
template<class T> struct vec_step<__vector_type<T, 4> > { enum { value = 4 }; };
template<class T> struct vec_step<__vector_type<T, 8> > { enum { value = 8 }; };
template<class T> struct vec_step<__vector_type<T, 16> > { enum { value = 16 }; };

#define vec_step(X)	(vec_step<typeof(X)>::value)

template<int N> struct __mask;
template<> struct __mask<2> { enum { value = 0x1U }; };
template<> struct __mask<3> { enum { value = 0x3U }; };
template<> struct __mask<4> { enum { value = 0x3U }; };
template<> struct __mask<8> { enum { value = 0x7U }; };
template<> struct __mask<16> { enum { value = 0xFU }; };

// shuffle
template<typename gentypem, typename ugentypen>
static inline typename __vector_type<typename __vector<gentypem>::basic_type, __vector<ugentypen>::components>::type
		shuffle(const gentypem &x, const ugentypen &mask)
{
	typename __vector_type<typename __vector<gentypem>::basic_type, __vector<ugentypen>::components>::type ret;
	const __size_t MASK = __mask<__vector<gentypem>::components>::value;
	for(__size_t i = 0 ; i < __vector<ugentypen>::components ; ++i)
		ret.v[i] = x.v[mask.v[i] & MASK];
	return ret;
}

// shuffle2
template<typename gentypem, typename ugentypen>
static inline typename __vector_type<typename __vector<gentypem>::basic_type, __vector<ugentypen>::components>::type
		shuffle2(const gentypem &x, const gentypem &y, const ugentypen &mask)
{
	typename __vector_type<typename __vector<gentypem>::basic_type, __vector<ugentypen>::components>::type ret;
	const __size_t MASK = __mask<__vector<gentypem>::components>::value;
	const __size_t MASK2 = MASK + 1U;
	for(__size_t i = 0 ; i < __vector<ugentypen>::components ; ++i)
		ret.v[i] = (mask.v[i] & MASK2) ? y.v[mask.v[i] & MASK] : x.v[mask.v[i] & MASK];
	return ret;
}

#endif

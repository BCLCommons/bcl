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
#ifndef __FREEOCL_OPENCL_C_PREINCLUDE_PRINTF_H__
#define __FREEOCL_OPENCL_C_PREINCLUDE_PRINTF_H__

#include <unistd.h>
#include <cstdarg>
#include "atomic.h"

#define FLAG_MINUS			0x01		// -
#define FLAG_PLUS			0x02		// +
#define FLAG_SPACE			0x04		// space
#define FLAG_ALTERNATIVE	0x08		// #
#define FLAG_0				0x10		// 0

#define LENGTH_L			0b00
#define LENGTH_H			0b01
#define LENGTH_HH			0b11
#define LENGTH_HL			0b10

static inline int __flag(int c)
{
	switch(c)
	{
	case '-':	return FLAG_MINUS;
	case '+':	return FLAG_PLUS;
	case ' ':	return FLAG_SPACE;
	case '#':	return FLAG_ALTERNATIVE;
	case '0':	return FLAG_0;
	default:
		return 0;
	}
}

static inline void __print_char(char c)
{
	write(STDOUT_FILENO, &c, 1);
}

static void __print_uint(unsigned long i, int base, int field_width, char padding_char, bool b_upper_case)
{
	char buf[20];
	int l = 0;
	if (i == 0)
	{
		for(int j = 0 ; j < field_width - 1 ; ++j)
			__print_char(padding_char);
		if (field_width != 0)
			__print_char('0');
		return;
	}
	const char digit_upper[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	const char digit_lower[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
	const char *digit = b_upper_case ? digit_upper : digit_lower;
	for( ; i != 0 ; i /= base, ++l)
		buf[l] = digit[i % base];
	for(int j = l ; j < field_width ; ++j)
		__print_char(padding_char);
	for(--l ; l >= 0 ; --l)
		__print_char(buf[l]);
}

static void __print_int(long i, int base, int field_width, char padding_char, bool b_upper_case)
{
	if (i < 0)
	{
		__print_char('-');
		i = -i;
	}
	__print_uint(i, base, field_width, padding_char, b_upper_case);
}

static void __print_double(double d, int precision, double base, bool b_upper_case)
{
	if (d < 0.0)
	{
		__print_char('-');
		d = -d;
	}
	if (isinf(d))
	{
		if (b_upper_case)
		{
			__print_char('I');
			__print_char('N');
			__print_char('F');
		}
		else
		{
			__print_char('i');
			__print_char('n');
			__print_char('f');
		}
		return;
	}
	if (isnan(d))
	{
		if (b_upper_case)
		{
			__print_char('N');
			__print_char('A');
			__print_char('N');
		}
		else
		{
			__print_char('n');
			__print_char('a');
			__print_char('n');
		}
		return;
	}
	precision = min(16, max(0, precision));
	double p = 1.0;
	for(int i = 0 ; i < precision ; ++i)
		p *= base;
	__print_uint(d, base, 1, '0', b_upper_case);
	__print_char('.');
	__print_uint((d - (unsigned long)d) * p + 0.5, base, precision, '0', b_upper_case);
}

static void __print_double_e(double d, int precision, bool b_upper_case)
{
	if (d < 0.0)
	{
		__print_char('-');
		d = -d;
	}
	if (isinf(d))
	{
		if (b_upper_case)
		{
			__print_char('I');
			__print_char('N');
			__print_char('F');
		}
		else
		{
			__print_char('i');
			__print_char('n');
			__print_char('f');
		}
		return;
	}
	if (isnan(d))
	{
		if (b_upper_case)
		{
			__print_char('N');
			__print_char('A');
			__print_char('N');
		}
		else
		{
			__print_char('n');
			__print_char('a');
			__print_char('n');
		}
		return;
	}
	if (d == 0.0)
	{
		__print_char('0');
		__print_char('.');
		for(int i = 0 ; i < precision ; ++i)
			__print_char('0');
		__print_char(b_upper_case ? 'E' : 'e');
		__print_char('+');
		__print_char('0');
		return;
	}
	int e = 0;
	double f = d;
	for(; f >= 10.0 ; f *= 0.1)
		e++;
	for(; f < 1.0 ; f *= 10.0)
		e--;
	__print_double(f, precision, 10.0, b_upper_case);
	__print_char(b_upper_case ? 'E' : 'e');
	if (e >= 0)
		__print_char('+');
	__print_int(e, 10, 1, '0', b_upper_case);
}

static void __print_double_a(double d, int precision, bool b_upper_case)
{
	if (d < 0.0)
	{
		__print_char('-');
		d = -d;
	}
	if (isinf(d))
	{
		if (b_upper_case)
		{
			__print_char('I');
			__print_char('N');
			__print_char('F');
		}
		else
		{
			__print_char('i');
			__print_char('n');
			__print_char('f');
		}
		return;
	}
	if (isnan(d))
	{
		if (b_upper_case)
		{
			__print_char('N');
			__print_char('A');
			__print_char('N');
		}
		else
		{
			__print_char('n');
			__print_char('a');
			__print_char('n');
		}
		return;
	}
	__print_char('0');
	__print_char(b_upper_case ? 'X' : 'x');
	if (d == 0.0)
	{
		__print_char('0');
		__print_char('.');
		for(int i = 0 ; i < precision ; ++i)
			__print_char('0');
		__print_char(b_upper_case ? 'P' : 'p');
		__print_char('+');
		__print_char('0');
		return;
	}
	int p = 0;
	double f = d;
	for(; f >= 16.0 ; f *= 0.5)
		p++;
	for(; f < 1.0 ; f *= 2.0)
		p--;
	__print_double(f, precision, 16.0, b_upper_case);
	__print_char(b_upper_case ? 'P' : 'p');
	if (p >= 0)
		__print_char('+');
	__print_int(p, 16, 1, '0', b_upper_case);
}

static void __print_double_g(double d, int precision, bool b_upper_case)
{
	if (isinf(d) || isnan(d) || d == 0.0)
	{
		__print_double(d, precision, 10.0, b_upper_case);
		return;
	}
	const int P = precision ? precision : 1;
	int e = 0;
	double f = d < 0.0 ? -d : d;
	for(; f >= 10.0 ; f *= 0.1)
		e++;
	for(; f < 1.0 ; f *= 10.0)
		e--;
	if (P > e && e >= -4)
		__print_double(d, P, 10.0, b_upper_case);
	else
		__print_double_e(d, P - 1, b_upper_case);
}

// Built-in printf function
static int printf(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	bool b_error = false;
	// Thread safety is not required by OpenCL 1.2 but it would make printf almost useless
	__atomic_lock();
	for(int i = 0 ; format[i] != 0 && !b_error ; ++i)
	{
		switch(format[i])
		{
		case '%':	// Directive
			{
				int flags = 0;
				for(++i ; __flag(format[i]) ; ++i)		// Parse flags
					flags |= __flag(format[i]);

				if (format[i] == 0)
				{
					b_error = true;
					continue;
				}

				int field_width = -1;
				if (format[i] > '0' && format[i] <= '9')
				{
					field_width = 0;
					for( ; format[i] != 0 && format[i] >= '0' && format[i] <= '9' ; ++i)		// Parse field width
						field_width = field_width * 10 + format[i] - '0';
				}

				int precision = -1;
				if (format[i] == '.')
				{
					precision = 0;
					for(++i ; format[i] != 0 && format[i] >= '0' && format[i] <= '9' ; ++i)		// Parse precision
						precision = precision * 10 + format[i] - '0';
				}

				int dimension = 1;
				const char vector_field_separator = ',';
				if (format[i] == 'v')
				{
					dimension = 0;
					for(++i ; format[i] != 0 && format[i] >= '0' && format[i] <= '9' ; ++i)	// Parse vector specifier
						dimension = dimension * 10 + format[i] - '0';
				}

				int length_modifier = -1;
				if (format[i] == 'h' || format[i] == 'l')			// Parse length modifier
				{
					length_modifier = (format[i] == 'h') ? LENGTH_H : LENGTH_L;
					++i;
					if (format[i] == 'h' || format[i] == 'l')
					{
						++i;
						if (length_modifier == LENGTH_L && format[i] == 'l')
						{
							b_error = true;
							continue;
						}
						length_modifier = (length_modifier << 1);
						length_modifier |= (format[i] == 'h') ? LENGTH_H : LENGTH_L;
					}
				}

				if (length_modifier == -1 && dimension > 1)
				{
					b_error = true;
					continue;
				}

#define READ_ARG(T)\
	switch(dimension)\
	{\
		case 2:\
			{\
				typedef T##2 _type;\
				const _type a = va_arg(ap, _type);\
				for(int j = 0 ; j < dimension ; ++j)\
				v[j] = a.v[j];\
			}\
			break;\
		case 3:\
			{\
				typedef T##3 _type;\
				const _type a = va_arg(ap, _type);\
				for(int j = 0 ; j < dimension ; ++j)\
					v[j] = a.v[j];\
			}\
			break;\
		case 4:\
			{\
				typedef T##4 _type;\
				const _type a = va_arg(ap, _type);\
				for(int j = 0 ; j < dimension ; ++j)\
					v[j] = a.v[j];\
			}\
			break;\
		case 8:\
			{\
				typedef T##8 _type;\
				const _type a = va_arg(ap, _type);\
				for(int j = 0 ; j < dimension ; ++j)\
					v[j] = a.v[j];\
			}\
			break;\
		case 16:\
			{\
				typedef T##16 _type;\
				const _type a = va_arg(ap, _type);\
				for(int j = 0 ; j < dimension ; ++j)\
					v[j] = a.v[j];\
			}\
			break;\
		}

				const int convertion_specifier = format[i];
				switch(convertion_specifier)
				{
				case 'd':
				case 'i':
					{
						long v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, long);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// char
								READ_ARG(__char)
								break;
							case LENGTH_H:		// short
								READ_ARG(__short)
								break;
							case LENGTH_HL:		// int
								READ_ARG(__int)
								break;
							case LENGTH_L:		// long
								READ_ARG(__long)
								break;
							}
						}
						char padding_char = '0';
						if (precision == -1)
						{
							if (field_width == -1)
								precision = 1;
							else
							{
								precision = field_width;
								padding_char = (flags & FLAG_0) ? '0' : ' ';
							}
						}
						for(int j = 0 ; j < dimension ; ++j)
						{
							if (v[j] >= 0)
							{
								if (flags & FLAG_PLUS)
									__print_char('+');
								else if (flags & FLAG_SPACE)
									__print_char(' ');
							}
							__print_int(v[j], 10, precision, padding_char, false);
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'o':
				case 'u':
				case 'x':
				case 'X':
					{
						unsigned long v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, unsigned long);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// uchar
								READ_ARG(__uchar)
								break;
							case LENGTH_H:		// ushort
								READ_ARG(__ushort)
								break;
							case LENGTH_HL:		// uint
								READ_ARG(__uint)
								break;
							case LENGTH_L:		// ulong
								READ_ARG(__ulong)
								break;
							}
						}
						const int base = convertion_specifier == 'o' ? 8
																	 : convertion_specifier == 'u' ? 10
																								   : 16;
						char padding_char = '0';
						if (precision == -1)
						{
							if (field_width == -1)
								precision = 1;
							else
							{
								precision = field_width;
								padding_char = (flags & FLAG_0) ? '0' : ' ';
							}
						}
						for(int j = 0 ; j < dimension ; ++j)
						{
							__print_uint(v[j], base, precision, padding_char, convertion_specifier == 'X');
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'f':
				case 'F':
					{
						double v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, double);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// 8bit float oO ?
							case LENGTH_H:		// halfn isn't supported either
								b_error = true;
								continue;
							case LENGTH_HL:		// float
								READ_ARG(__float)
								break;
							case LENGTH_L:		// double
								READ_ARG(__double)
								break;
							}
						}
						if (precision == -1)
							precision = 6;
						for(int j = 0 ; j < dimension ; ++j)
						{
							if (v[j] >= 0.0)
							{
								if (flags & FLAG_PLUS)
									__print_char('+');
								else if (flags & FLAG_SPACE)
									__print_char(' ');
							}
							__print_double(v[j], precision, 10.0, convertion_specifier == 'F');
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'e':
				case 'E':
					{
						double v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, double);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// 8bit float oO ?
							case LENGTH_H:		// halfn isn't supported either
								b_error = true;
								continue;
							case LENGTH_HL:		// float
								READ_ARG(__float)
								break;
							case LENGTH_L:		// double
								READ_ARG(__double)
								break;
							}
						}
						if (precision == -1)
							precision = 6;
						for(int j = 0 ; j < dimension ; ++j)
						{
							if (v[j] >= 0.0)
							{
								if (flags & FLAG_PLUS)
									__print_char('+');
								else if (flags & FLAG_SPACE)
									__print_char(' ');
							}
							__print_double_e(v[j], precision, convertion_specifier == 'E');
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'g':
				case 'G':
					{
						double v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, double);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// 8bit float oO ?
							case LENGTH_H:		// halfn isn't supported either
								b_error = true;
								continue;
							case LENGTH_HL:		// float
								READ_ARG(__float)
								break;
							case LENGTH_L:		// double
								READ_ARG(__double)
								break;
							}
						}
						if (precision == -1)
							precision = 6;
						for(int j = 0 ; j < dimension ; ++j)
						{
							if (v[j] >= 0.0)
							{
								if (flags & FLAG_PLUS)
									__print_char('+');
								else if (flags & FLAG_SPACE)
									__print_char(' ');
							}
							__print_double_g(v[j], precision, convertion_specifier == 'G');
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'a':
				case 'A':
					{
						double v[16];
						if (dimension == 1)
							v[0] = va_arg(ap, double);
						else
						{
							switch(length_modifier)
							{
							case LENGTH_HH:		// 8bit float oO ?
							case LENGTH_H:		// halfn isn't supported either
								b_error = true;
								continue;
							case LENGTH_HL:		// float
								READ_ARG(__float)
								break;
							case LENGTH_L:		// double
								READ_ARG(__double)
								break;
							}
						}
						if (precision == -1)
							precision = 6;
						for(int j = 0 ; j < dimension ; ++j)
						{
							if (v[j] >= 0.0)
							{
								if (flags & FLAG_PLUS)
									__print_char('+');
								else if (flags & FLAG_SPACE)
									__print_char(' ');
							}
							__print_double_a(v[j], precision, convertion_specifier == 'A');
							if (j + 1 < dimension)
								__print_char(vector_field_separator);
						}
					}
					break;
				case 'c':
					__print_char(va_arg(ap, int));
					break;
				case 's':
					{
						const char *s = va_arg(ap, const char*);
						for(size_t j = 0 ; s[j] != 0 ; ++j)
							__print_char(s[j]);
					}
					break;
				case 'p':
					__print_char('0');
					__print_char('x');
					__print_uint((__size_t)va_arg(ap, const void*), 16, sizeof(void*) * 2, '0', false);
					break;
				case '%':
					__print_char('%');
					break;
				default:
					b_error = true;
					continue;
				}
			}
			break;
		default:
			__print_char(format[i]);
		}
	}

	__atomic_unlock();

	va_end(ap);
	return b_error ? -1 : 0;
}

#endif

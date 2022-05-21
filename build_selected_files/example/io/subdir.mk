################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/io/example_io.cpp \
../example/io/example_io_binary_serialize.cpp \
../example/io/example_io_directory.cpp \
../example/io/example_io_directory_entry.cpp \
../example/io/example_io_file.cpp \
../example/io/example_io_file_stream_buffer.cpp \
../example/io/example_io_fixed_line_width_writer.cpp \
../example/io/example_io_ifstream.cpp \
../example/io/example_io_ofstream.cpp \
../example/io/example_io_serialization_builtin.cpp \
../example/io/example_io_serialization_container.cpp \
../example/io/example_io_serialization_map.cpp \
../example/io/example_io_serialization_with_check.cpp \
../example/io/example_io_serialization_with_min_max.cpp \
../example/io/example_io_serialization_wrapper.cpp \
../example/io/example_io_serialize.cpp \
../example/io/example_io_stream_buffer_classes.cpp 

OBJS += \
./example/io/example_io.o \
./example/io/example_io_binary_serialize.o \
./example/io/example_io_directory.o \
./example/io/example_io_directory_entry.o \
./example/io/example_io_file.o \
./example/io/example_io_file_stream_buffer.o \
./example/io/example_io_fixed_line_width_writer.o \
./example/io/example_io_ifstream.o \
./example/io/example_io_ofstream.o \
./example/io/example_io_serialization_builtin.o \
./example/io/example_io_serialization_container.o \
./example/io/example_io_serialization_map.o \
./example/io/example_io_serialization_with_check.o \
./example/io/example_io_serialization_with_min_max.o \
./example/io/example_io_serialization_wrapper.o \
./example/io/example_io_serialize.o \
./example/io/example_io_stream_buffer_classes.o 

CPP_DEPS += \
./example/io/example_io.d \
./example/io/example_io_binary_serialize.d \
./example/io/example_io_directory.d \
./example/io/example_io_directory_entry.d \
./example/io/example_io_file.d \
./example/io/example_io_file_stream_buffer.d \
./example/io/example_io_fixed_line_width_writer.d \
./example/io/example_io_ifstream.d \
./example/io/example_io_ofstream.d \
./example/io/example_io_serialization_builtin.d \
./example/io/example_io_serialization_container.d \
./example/io/example_io_serialization_map.d \
./example/io/example_io_serialization_with_check.d \
./example/io/example_io_serialization_with_min_max.d \
./example/io/example_io_serialization_wrapper.d \
./example/io/example_io_serialize.d \
./example/io/example_io_stream_buffer_classes.d 


# Each subdirectory must supply rules for building sources it contributes
example/io/%.o: ../example/io/%.cpp example/io/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/io/example_io.cpp \
../bcl_test_release/bcl/example/io/example_io_binary_serialize.cpp \
../bcl_test_release/bcl/example/io/example_io_directory.cpp \
../bcl_test_release/bcl/example/io/example_io_directory_entry.cpp \
../bcl_test_release/bcl/example/io/example_io_file.cpp \
../bcl_test_release/bcl/example/io/example_io_file_stream_buffer.cpp \
../bcl_test_release/bcl/example/io/example_io_fixed_line_width_writer.cpp \
../bcl_test_release/bcl/example/io/example_io_ifstream.cpp \
../bcl_test_release/bcl/example/io/example_io_ofstream.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_builtin.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_container.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_map.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_with_check.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_with_min_max.cpp \
../bcl_test_release/bcl/example/io/example_io_serialization_wrapper.cpp \
../bcl_test_release/bcl/example/io/example_io_serialize.cpp \
../bcl_test_release/bcl/example/io/example_io_stream_buffer_classes.cpp 

OBJS += \
./bcl_test_release/bcl/example/io/example_io.o \
./bcl_test_release/bcl/example/io/example_io_binary_serialize.o \
./bcl_test_release/bcl/example/io/example_io_directory.o \
./bcl_test_release/bcl/example/io/example_io_directory_entry.o \
./bcl_test_release/bcl/example/io/example_io_file.o \
./bcl_test_release/bcl/example/io/example_io_file_stream_buffer.o \
./bcl_test_release/bcl/example/io/example_io_fixed_line_width_writer.o \
./bcl_test_release/bcl/example/io/example_io_ifstream.o \
./bcl_test_release/bcl/example/io/example_io_ofstream.o \
./bcl_test_release/bcl/example/io/example_io_serialization_builtin.o \
./bcl_test_release/bcl/example/io/example_io_serialization_container.o \
./bcl_test_release/bcl/example/io/example_io_serialization_map.o \
./bcl_test_release/bcl/example/io/example_io_serialization_with_check.o \
./bcl_test_release/bcl/example/io/example_io_serialization_with_min_max.o \
./bcl_test_release/bcl/example/io/example_io_serialization_wrapper.o \
./bcl_test_release/bcl/example/io/example_io_serialize.o \
./bcl_test_release/bcl/example/io/example_io_stream_buffer_classes.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/io/example_io.d \
./bcl_test_release/bcl/example/io/example_io_binary_serialize.d \
./bcl_test_release/bcl/example/io/example_io_directory.d \
./bcl_test_release/bcl/example/io/example_io_directory_entry.d \
./bcl_test_release/bcl/example/io/example_io_file.d \
./bcl_test_release/bcl/example/io/example_io_file_stream_buffer.d \
./bcl_test_release/bcl/example/io/example_io_fixed_line_width_writer.d \
./bcl_test_release/bcl/example/io/example_io_ifstream.d \
./bcl_test_release/bcl/example/io/example_io_ofstream.d \
./bcl_test_release/bcl/example/io/example_io_serialization_builtin.d \
./bcl_test_release/bcl/example/io/example_io_serialization_container.d \
./bcl_test_release/bcl/example/io/example_io_serialization_map.d \
./bcl_test_release/bcl/example/io/example_io_serialization_with_check.d \
./bcl_test_release/bcl/example/io/example_io_serialization_with_min_max.d \
./bcl_test_release/bcl/example/io/example_io_serialization_wrapper.d \
./bcl_test_release/bcl/example/io/example_io_serialize.d \
./bcl_test_release/bcl/example/io/example_io_stream_buffer_classes.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/io/%.o: ../bcl_test_release/bcl/example/io/%.cpp bcl_test_release/bcl/example/io/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



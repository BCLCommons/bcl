################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/io/bcl_io.cpp \
../source/io/bcl_io_binary_serialize.cpp \
../source/io/bcl_io_directory.cpp \
../source/io/bcl_io_directory_entry.cpp \
../source/io/bcl_io_file.cpp \
../source/io/bcl_io_file_stream_buffer.cpp \
../source/io/bcl_io_file_stream_buffer_encrypted.cpp \
../source/io/bcl_io_fixed_line_width_writer.cpp \
../source/io/bcl_io_ifstream.cpp \
../source/io/bcl_io_ofstream.cpp \
../source/io/bcl_io_serialization.cpp \
../source/io/bcl_io_serialization_builtin.cpp \
../source/io/bcl_io_serialization_interface.cpp \
../source/io/bcl_io_serialization_with_check.cpp \
../source/io/bcl_io_serialization_with_min_max.cpp \
../source/io/bcl_io_serialize.cpp \
../source/io/bcl_io_serializer.cpp \
../source/io/bcl_io_stream_buffer_classes.cpp \
../source/io/bcl_io_stream_buffer_interface.cpp \
../source/io/bcl_io_stream_interface.cpp \
../source/io/bcl_io_validation_result.cpp 

OBJS += \
./source/io/bcl_io.o \
./source/io/bcl_io_binary_serialize.o \
./source/io/bcl_io_directory.o \
./source/io/bcl_io_directory_entry.o \
./source/io/bcl_io_file.o \
./source/io/bcl_io_file_stream_buffer.o \
./source/io/bcl_io_file_stream_buffer_encrypted.o \
./source/io/bcl_io_fixed_line_width_writer.o \
./source/io/bcl_io_ifstream.o \
./source/io/bcl_io_ofstream.o \
./source/io/bcl_io_serialization.o \
./source/io/bcl_io_serialization_builtin.o \
./source/io/bcl_io_serialization_interface.o \
./source/io/bcl_io_serialization_with_check.o \
./source/io/bcl_io_serialization_with_min_max.o \
./source/io/bcl_io_serialize.o \
./source/io/bcl_io_serializer.o \
./source/io/bcl_io_stream_buffer_classes.o \
./source/io/bcl_io_stream_buffer_interface.o \
./source/io/bcl_io_stream_interface.o \
./source/io/bcl_io_validation_result.o 

CPP_DEPS += \
./source/io/bcl_io.d \
./source/io/bcl_io_binary_serialize.d \
./source/io/bcl_io_directory.d \
./source/io/bcl_io_directory_entry.d \
./source/io/bcl_io_file.d \
./source/io/bcl_io_file_stream_buffer.d \
./source/io/bcl_io_file_stream_buffer_encrypted.d \
./source/io/bcl_io_fixed_line_width_writer.d \
./source/io/bcl_io_ifstream.d \
./source/io/bcl_io_ofstream.d \
./source/io/bcl_io_serialization.d \
./source/io/bcl_io_serialization_builtin.d \
./source/io/bcl_io_serialization_interface.d \
./source/io/bcl_io_serialization_with_check.d \
./source/io/bcl_io_serialization_with_min_max.d \
./source/io/bcl_io_serialize.d \
./source/io/bcl_io_serializer.d \
./source/io/bcl_io_stream_buffer_classes.d \
./source/io/bcl_io_stream_buffer_interface.d \
./source/io/bcl_io_stream_interface.d \
./source/io/bcl_io_validation_result.d 


# Each subdirectory must supply rules for building sources it contributes
source/io/%.o: ../source/io/%.cpp source/io/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



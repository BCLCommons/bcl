################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/linal/bcl_linal.cpp \
../source/linal/bcl_linal_distance_geometry.cpp \
../source/linal/bcl_linal_matrix.cpp \
../source/linal/bcl_linal_matrix2x2.cpp \
../source/linal/bcl_linal_matrix3x3.cpp \
../source/linal/bcl_linal_matrix_const_interface.cpp \
../source/linal/bcl_linal_matrix_interface.cpp \
../source/linal/bcl_linal_matrix_inversion_cholesky.cpp \
../source/linal/bcl_linal_matrix_inversion_gauss_jordan.cpp \
../source/linal/bcl_linal_matrix_inversion_interface.cpp \
../source/linal/bcl_linal_matrix_inversion_moore_penrose.cpp \
../source/linal/bcl_linal_operations.cpp \
../source/linal/bcl_linal_operations_cpu.cpp \
../source/linal/bcl_linal_operations_interface.cpp \
../source/linal/bcl_linal_vector.cpp \
../source/linal/bcl_linal_vector_2d.cpp \
../source/linal/bcl_linal_vector_2d_operations.cpp \
../source/linal/bcl_linal_vector_3d.cpp \
../source/linal/bcl_linal_vector_3d_operations.cpp \
../source/linal/bcl_linal_vector_const_interface.cpp \
../source/linal/bcl_linal_vector_interface.cpp 

OBJS += \
./source/linal/bcl_linal.o \
./source/linal/bcl_linal_distance_geometry.o \
./source/linal/bcl_linal_matrix.o \
./source/linal/bcl_linal_matrix2x2.o \
./source/linal/bcl_linal_matrix3x3.o \
./source/linal/bcl_linal_matrix_const_interface.o \
./source/linal/bcl_linal_matrix_interface.o \
./source/linal/bcl_linal_matrix_inversion_cholesky.o \
./source/linal/bcl_linal_matrix_inversion_gauss_jordan.o \
./source/linal/bcl_linal_matrix_inversion_interface.o \
./source/linal/bcl_linal_matrix_inversion_moore_penrose.o \
./source/linal/bcl_linal_operations.o \
./source/linal/bcl_linal_operations_cpu.o \
./source/linal/bcl_linal_operations_interface.o \
./source/linal/bcl_linal_vector.o \
./source/linal/bcl_linal_vector_2d.o \
./source/linal/bcl_linal_vector_2d_operations.o \
./source/linal/bcl_linal_vector_3d.o \
./source/linal/bcl_linal_vector_3d_operations.o \
./source/linal/bcl_linal_vector_const_interface.o \
./source/linal/bcl_linal_vector_interface.o 

CPP_DEPS += \
./source/linal/bcl_linal.d \
./source/linal/bcl_linal_distance_geometry.d \
./source/linal/bcl_linal_matrix.d \
./source/linal/bcl_linal_matrix2x2.d \
./source/linal/bcl_linal_matrix3x3.d \
./source/linal/bcl_linal_matrix_const_interface.d \
./source/linal/bcl_linal_matrix_interface.d \
./source/linal/bcl_linal_matrix_inversion_cholesky.d \
./source/linal/bcl_linal_matrix_inversion_gauss_jordan.d \
./source/linal/bcl_linal_matrix_inversion_interface.d \
./source/linal/bcl_linal_matrix_inversion_moore_penrose.d \
./source/linal/bcl_linal_operations.d \
./source/linal/bcl_linal_operations_cpu.d \
./source/linal/bcl_linal_operations_interface.d \
./source/linal/bcl_linal_vector.d \
./source/linal/bcl_linal_vector_2d.d \
./source/linal/bcl_linal_vector_2d_operations.d \
./source/linal/bcl_linal_vector_3d.d \
./source/linal/bcl_linal_vector_3d_operations.d \
./source/linal/bcl_linal_vector_const_interface.d \
./source/linal/bcl_linal_vector_interface.d 


# Each subdirectory must supply rules for building sources it contributes
source/linal/%.o: ../source/linal/%.cpp source/linal/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/linal/example_linal.cpp \
../example/linal/example_linal_distance_geometry.cpp \
../example/linal/example_linal_matrix.cpp \
../example/linal/example_linal_matrix3x3.cpp \
../example/linal/example_linal_matrix_const_reference.cpp \
../example/linal/example_linal_matrix_inversion_cholesky.cpp \
../example/linal/example_linal_matrix_inversion_gauss_jordan.cpp \
../example/linal/example_linal_matrix_inversion_interface.cpp \
../example/linal/example_linal_matrix_inversion_moore_penrose.cpp \
../example/linal/example_linal_matrix_operations.cpp \
../example/linal/example_linal_matrix_reference.cpp \
../example/linal/example_linal_operations.cpp \
../example/linal/example_linal_principal_component_analysis.cpp \
../example/linal/example_linal_symmetric_eigensolver.cpp \
../example/linal/example_linal_vector.cpp \
../example/linal/example_linal_vector_2d.cpp \
../example/linal/example_linal_vector_2d_operations.cpp \
../example/linal/example_linal_vector_3d.cpp \
../example/linal/example_linal_vector_3d_operations.cpp \
../example/linal/example_linal_vector_const_reference.cpp \
../example/linal/example_linal_vector_nd.cpp \
../example/linal/example_linal_vector_operations.cpp \
../example/linal/example_linal_vector_reference.cpp 

OBJS += \
./example/linal/example_linal.o \
./example/linal/example_linal_distance_geometry.o \
./example/linal/example_linal_matrix.o \
./example/linal/example_linal_matrix3x3.o \
./example/linal/example_linal_matrix_const_reference.o \
./example/linal/example_linal_matrix_inversion_cholesky.o \
./example/linal/example_linal_matrix_inversion_gauss_jordan.o \
./example/linal/example_linal_matrix_inversion_interface.o \
./example/linal/example_linal_matrix_inversion_moore_penrose.o \
./example/linal/example_linal_matrix_operations.o \
./example/linal/example_linal_matrix_reference.o \
./example/linal/example_linal_operations.o \
./example/linal/example_linal_principal_component_analysis.o \
./example/linal/example_linal_symmetric_eigensolver.o \
./example/linal/example_linal_vector.o \
./example/linal/example_linal_vector_2d.o \
./example/linal/example_linal_vector_2d_operations.o \
./example/linal/example_linal_vector_3d.o \
./example/linal/example_linal_vector_3d_operations.o \
./example/linal/example_linal_vector_const_reference.o \
./example/linal/example_linal_vector_nd.o \
./example/linal/example_linal_vector_operations.o \
./example/linal/example_linal_vector_reference.o 

CPP_DEPS += \
./example/linal/example_linal.d \
./example/linal/example_linal_distance_geometry.d \
./example/linal/example_linal_matrix.d \
./example/linal/example_linal_matrix3x3.d \
./example/linal/example_linal_matrix_const_reference.d \
./example/linal/example_linal_matrix_inversion_cholesky.d \
./example/linal/example_linal_matrix_inversion_gauss_jordan.d \
./example/linal/example_linal_matrix_inversion_interface.d \
./example/linal/example_linal_matrix_inversion_moore_penrose.d \
./example/linal/example_linal_matrix_operations.d \
./example/linal/example_linal_matrix_reference.d \
./example/linal/example_linal_operations.d \
./example/linal/example_linal_principal_component_analysis.d \
./example/linal/example_linal_symmetric_eigensolver.d \
./example/linal/example_linal_vector.d \
./example/linal/example_linal_vector_2d.d \
./example/linal/example_linal_vector_2d_operations.d \
./example/linal/example_linal_vector_3d.d \
./example/linal/example_linal_vector_3d_operations.d \
./example/linal/example_linal_vector_const_reference.d \
./example/linal/example_linal_vector_nd.d \
./example/linal/example_linal_vector_operations.d \
./example/linal/example_linal_vector_reference.d 


# Each subdirectory must supply rules for building sources it contributes
example/linal/%.o: ../example/linal/%.cpp example/linal/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



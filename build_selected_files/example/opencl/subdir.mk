################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/opencl/example_opencl.cpp \
../example/opencl/example_opencl_approximator_resilient_propagation.cpp \
../example/opencl/example_opencl_approximator_sequential_minimial_optimization.cpp \
../example/opencl/example_opencl_approximator_simple_propagation.cpp \
../example/opencl/example_opencl_arg_max.cpp \
../example/opencl/example_opencl_arg_min.cpp \
../example/opencl/example_opencl_dataset_min_max.cpp \
../example/opencl/example_opencl_device.cpp \
../example/opencl/example_opencl_euclidean_distance.cpp \
../example/opencl/example_opencl_extension_data.cpp \
../example/opencl/example_opencl_feature_similarity_measures.cpp \
../example/opencl/example_opencl_insertion_sort.cpp \
../example/opencl/example_opencl_kappa_nearest_neighbor.cpp \
../example/opencl/example_opencl_matrix.cpp \
../example/opencl/example_opencl_matrix_add.cpp \
../example/opencl/example_opencl_matrix_multiply.cpp \
../example/opencl/example_opencl_matrix_transpose.cpp \
../example/opencl/example_opencl_neural_network.cpp \
../example/opencl/example_opencl_platform.cpp \
../example/opencl/example_opencl_quality_gdt.cpp \
../example/opencl/example_opencl_quality_lcs.cpp \
../example/opencl/example_opencl_quality_rmsd.cpp \
../example/opencl/example_opencl_rmsd.cpp \
../example/opencl/example_opencl_saxs_debye.cpp \
../example/opencl/example_opencl_singular_value_decomposition.cpp \
../example/opencl/example_opencl_support_vector_machine.cpp \
../example/opencl/example_opencl_tools.cpp \
../example/opencl/example_opencl_transfer_function_gaussian.cpp \
../example/opencl/example_opencl_transfer_function_sigmoid.cpp \
../example/opencl/example_opencl_vector.cpp \
../example/opencl/example_opencl_vector_matrix_add.cpp 

OBJS += \
./example/opencl/example_opencl.o \
./example/opencl/example_opencl_approximator_resilient_propagation.o \
./example/opencl/example_opencl_approximator_sequential_minimial_optimization.o \
./example/opencl/example_opencl_approximator_simple_propagation.o \
./example/opencl/example_opencl_arg_max.o \
./example/opencl/example_opencl_arg_min.o \
./example/opencl/example_opencl_dataset_min_max.o \
./example/opencl/example_opencl_device.o \
./example/opencl/example_opencl_euclidean_distance.o \
./example/opencl/example_opencl_extension_data.o \
./example/opencl/example_opencl_feature_similarity_measures.o \
./example/opencl/example_opencl_insertion_sort.o \
./example/opencl/example_opencl_kappa_nearest_neighbor.o \
./example/opencl/example_opencl_matrix.o \
./example/opencl/example_opencl_matrix_add.o \
./example/opencl/example_opencl_matrix_multiply.o \
./example/opencl/example_opencl_matrix_transpose.o \
./example/opencl/example_opencl_neural_network.o \
./example/opencl/example_opencl_platform.o \
./example/opencl/example_opencl_quality_gdt.o \
./example/opencl/example_opencl_quality_lcs.o \
./example/opencl/example_opencl_quality_rmsd.o \
./example/opencl/example_opencl_rmsd.o \
./example/opencl/example_opencl_saxs_debye.o \
./example/opencl/example_opencl_singular_value_decomposition.o \
./example/opencl/example_opencl_support_vector_machine.o \
./example/opencl/example_opencl_tools.o \
./example/opencl/example_opencl_transfer_function_gaussian.o \
./example/opencl/example_opencl_transfer_function_sigmoid.o \
./example/opencl/example_opencl_vector.o \
./example/opencl/example_opencl_vector_matrix_add.o 

CPP_DEPS += \
./example/opencl/example_opencl.d \
./example/opencl/example_opencl_approximator_resilient_propagation.d \
./example/opencl/example_opencl_approximator_sequential_minimial_optimization.d \
./example/opencl/example_opencl_approximator_simple_propagation.d \
./example/opencl/example_opencl_arg_max.d \
./example/opencl/example_opencl_arg_min.d \
./example/opencl/example_opencl_dataset_min_max.d \
./example/opencl/example_opencl_device.d \
./example/opencl/example_opencl_euclidean_distance.d \
./example/opencl/example_opencl_extension_data.d \
./example/opencl/example_opencl_feature_similarity_measures.d \
./example/opencl/example_opencl_insertion_sort.d \
./example/opencl/example_opencl_kappa_nearest_neighbor.d \
./example/opencl/example_opencl_matrix.d \
./example/opencl/example_opencl_matrix_add.d \
./example/opencl/example_opencl_matrix_multiply.d \
./example/opencl/example_opencl_matrix_transpose.d \
./example/opencl/example_opencl_neural_network.d \
./example/opencl/example_opencl_platform.d \
./example/opencl/example_opencl_quality_gdt.d \
./example/opencl/example_opencl_quality_lcs.d \
./example/opencl/example_opencl_quality_rmsd.d \
./example/opencl/example_opencl_rmsd.d \
./example/opencl/example_opencl_saxs_debye.d \
./example/opencl/example_opencl_singular_value_decomposition.d \
./example/opencl/example_opencl_support_vector_machine.d \
./example/opencl/example_opencl_tools.d \
./example/opencl/example_opencl_transfer_function_gaussian.d \
./example/opencl/example_opencl_transfer_function_sigmoid.d \
./example/opencl/example_opencl_vector.d \
./example/opencl/example_opencl_vector_matrix_add.d 


# Each subdirectory must supply rules for building sources it contributes
example/opencl/%.o: ../example/opencl/%.cpp example/opencl/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



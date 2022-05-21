################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/opencl/bcl_opencl.cpp \
../source/opencl/bcl_opencl_approximator_resilient_propagation.cpp \
../source/opencl/bcl_opencl_approximator_sequential_minimial_optimization.cpp \
../source/opencl/bcl_opencl_approximator_simple_propagation.cpp \
../source/opencl/bcl_opencl_buffer.cpp \
../source/opencl/bcl_opencl_command_queue.cpp \
../source/opencl/bcl_opencl_context.cpp \
../source/opencl/bcl_opencl_density_fit_protein_minimizer_powell.cpp \
../source/opencl/bcl_opencl_density_simulate_gaussian_sphere.cpp \
../source/opencl/bcl_opencl_device.cpp \
../source/opencl/bcl_opencl_extension_data.cpp \
../source/opencl/bcl_opencl_extensions.cpp \
../source/opencl/bcl_opencl_feature_similarity_measures.cpp \
../source/opencl/bcl_opencl_insertion_sort.cpp \
../source/opencl/bcl_opencl_kappa_nearest_neighbor.cpp \
../source/opencl/bcl_opencl_kernel_source_alternative.cpp \
../source/opencl/bcl_opencl_kernel_source_file.cpp \
../source/opencl/bcl_opencl_kernel_source_interface.cpp \
../source/opencl/bcl_opencl_kernel_source_string.cpp \
../source/opencl/bcl_opencl_kernel_sources.cpp \
../source/opencl/bcl_opencl_matrix3x3.cpp \
../source/opencl/bcl_opencl_neural_network.cpp \
../source/opencl/bcl_opencl_operations.cpp \
../source/opencl/bcl_opencl_platform.cpp \
../source/opencl/bcl_opencl_protein_agreement_ccc.cpp \
../source/opencl/bcl_opencl_quality_gdt.cpp \
../source/opencl/bcl_opencl_quality_lcs.cpp \
../source/opencl/bcl_opencl_quality_rmsd.cpp \
../source/opencl/bcl_opencl_rmsd.cpp \
../source/opencl/bcl_opencl_saxs_debye.cpp \
../source/opencl/bcl_opencl_support_vector_machine.cpp \
../source/opencl/bcl_opencl_tools.cpp 

OBJS += \
./source/opencl/bcl_opencl.o \
./source/opencl/bcl_opencl_approximator_resilient_propagation.o \
./source/opencl/bcl_opencl_approximator_sequential_minimial_optimization.o \
./source/opencl/bcl_opencl_approximator_simple_propagation.o \
./source/opencl/bcl_opencl_buffer.o \
./source/opencl/bcl_opencl_command_queue.o \
./source/opencl/bcl_opencl_context.o \
./source/opencl/bcl_opencl_density_fit_protein_minimizer_powell.o \
./source/opencl/bcl_opencl_density_simulate_gaussian_sphere.o \
./source/opencl/bcl_opencl_device.o \
./source/opencl/bcl_opencl_extension_data.o \
./source/opencl/bcl_opencl_extensions.o \
./source/opencl/bcl_opencl_feature_similarity_measures.o \
./source/opencl/bcl_opencl_insertion_sort.o \
./source/opencl/bcl_opencl_kappa_nearest_neighbor.o \
./source/opencl/bcl_opencl_kernel_source_alternative.o \
./source/opencl/bcl_opencl_kernel_source_file.o \
./source/opencl/bcl_opencl_kernel_source_interface.o \
./source/opencl/bcl_opencl_kernel_source_string.o \
./source/opencl/bcl_opencl_kernel_sources.o \
./source/opencl/bcl_opencl_matrix3x3.o \
./source/opencl/bcl_opencl_neural_network.o \
./source/opencl/bcl_opencl_operations.o \
./source/opencl/bcl_opencl_platform.o \
./source/opencl/bcl_opencl_protein_agreement_ccc.o \
./source/opencl/bcl_opencl_quality_gdt.o \
./source/opencl/bcl_opencl_quality_lcs.o \
./source/opencl/bcl_opencl_quality_rmsd.o \
./source/opencl/bcl_opencl_rmsd.o \
./source/opencl/bcl_opencl_saxs_debye.o \
./source/opencl/bcl_opencl_support_vector_machine.o \
./source/opencl/bcl_opencl_tools.o 

CPP_DEPS += \
./source/opencl/bcl_opencl.d \
./source/opencl/bcl_opencl_approximator_resilient_propagation.d \
./source/opencl/bcl_opencl_approximator_sequential_minimial_optimization.d \
./source/opencl/bcl_opencl_approximator_simple_propagation.d \
./source/opencl/bcl_opencl_buffer.d \
./source/opencl/bcl_opencl_command_queue.d \
./source/opencl/bcl_opencl_context.d \
./source/opencl/bcl_opencl_density_fit_protein_minimizer_powell.d \
./source/opencl/bcl_opencl_density_simulate_gaussian_sphere.d \
./source/opencl/bcl_opencl_device.d \
./source/opencl/bcl_opencl_extension_data.d \
./source/opencl/bcl_opencl_extensions.d \
./source/opencl/bcl_opencl_feature_similarity_measures.d \
./source/opencl/bcl_opencl_insertion_sort.d \
./source/opencl/bcl_opencl_kappa_nearest_neighbor.d \
./source/opencl/bcl_opencl_kernel_source_alternative.d \
./source/opencl/bcl_opencl_kernel_source_file.d \
./source/opencl/bcl_opencl_kernel_source_interface.d \
./source/opencl/bcl_opencl_kernel_source_string.d \
./source/opencl/bcl_opencl_kernel_sources.d \
./source/opencl/bcl_opencl_matrix3x3.d \
./source/opencl/bcl_opencl_neural_network.d \
./source/opencl/bcl_opencl_operations.d \
./source/opencl/bcl_opencl_platform.d \
./source/opencl/bcl_opencl_protein_agreement_ccc.d \
./source/opencl/bcl_opencl_quality_gdt.d \
./source/opencl/bcl_opencl_quality_lcs.d \
./source/opencl/bcl_opencl_quality_rmsd.d \
./source/opencl/bcl_opencl_rmsd.d \
./source/opencl/bcl_opencl_saxs_debye.d \
./source/opencl/bcl_opencl_support_vector_machine.d \
./source/opencl/bcl_opencl_tools.d 


# Each subdirectory must supply rules for building sources it contributes
source/opencl/%.o: ../source/opencl/%.cpp source/opencl/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



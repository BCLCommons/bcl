################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/density/example_density.cpp \
../bcl_test_release/bcl/example/density/example_density_connectivity.cpp \
../bcl_test_release/bcl/example/density/example_density_map.cpp \
../bcl_test_release/bcl/example/density/example_density_map_cylindrical.cpp \
../bcl_test_release/bcl/example/density/example_density_mask_3d.cpp \
../bcl_test_release/bcl/example/density/example_density_protein_agreement_ccc.cpp \
../bcl_test_release/bcl/example/density/example_density_protein_agreement_likelihood.cpp \
../bcl_test_release/bcl/example/density/example_density_protein_agreements.cpp \
../bcl_test_release/bcl/example/density/example_density_simulate_default.cpp \
../bcl_test_release/bcl/example/density/example_density_simulate_gaussian_sphere.cpp \
../bcl_test_release/bcl/example/density/example_density_simulators.cpp 

OBJS += \
./bcl_test_release/bcl/example/density/example_density.o \
./bcl_test_release/bcl/example/density/example_density_connectivity.o \
./bcl_test_release/bcl/example/density/example_density_map.o \
./bcl_test_release/bcl/example/density/example_density_map_cylindrical.o \
./bcl_test_release/bcl/example/density/example_density_mask_3d.o \
./bcl_test_release/bcl/example/density/example_density_protein_agreement_ccc.o \
./bcl_test_release/bcl/example/density/example_density_protein_agreement_likelihood.o \
./bcl_test_release/bcl/example/density/example_density_protein_agreements.o \
./bcl_test_release/bcl/example/density/example_density_simulate_default.o \
./bcl_test_release/bcl/example/density/example_density_simulate_gaussian_sphere.o \
./bcl_test_release/bcl/example/density/example_density_simulators.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/density/example_density.d \
./bcl_test_release/bcl/example/density/example_density_connectivity.d \
./bcl_test_release/bcl/example/density/example_density_map.d \
./bcl_test_release/bcl/example/density/example_density_map_cylindrical.d \
./bcl_test_release/bcl/example/density/example_density_mask_3d.d \
./bcl_test_release/bcl/example/density/example_density_protein_agreement_ccc.d \
./bcl_test_release/bcl/example/density/example_density_protein_agreement_likelihood.d \
./bcl_test_release/bcl/example/density/example_density_protein_agreements.d \
./bcl_test_release/bcl/example/density/example_density_simulate_default.d \
./bcl_test_release/bcl/example/density/example_density_simulate_gaussian_sphere.d \
./bcl_test_release/bcl/example/density/example_density_simulators.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/density/%.o: ../bcl_test_release/bcl/example/density/%.cpp bcl_test_release/bcl/example/density/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



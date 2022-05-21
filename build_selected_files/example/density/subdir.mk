################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/density/example_density.cpp \
../example/density/example_density_connectivity.cpp \
../example/density/example_density_map.cpp \
../example/density/example_density_map_cylindrical.cpp \
../example/density/example_density_mask_3d.cpp \
../example/density/example_density_protein_agreement_ccc.cpp \
../example/density/example_density_protein_agreement_likelihood.cpp \
../example/density/example_density_protein_agreements.cpp \
../example/density/example_density_simulate_default.cpp \
../example/density/example_density_simulate_gaussian_sphere.cpp \
../example/density/example_density_simulators.cpp 

OBJS += \
./example/density/example_density.o \
./example/density/example_density_connectivity.o \
./example/density/example_density_map.o \
./example/density/example_density_map_cylindrical.o \
./example/density/example_density_mask_3d.o \
./example/density/example_density_protein_agreement_ccc.o \
./example/density/example_density_protein_agreement_likelihood.o \
./example/density/example_density_protein_agreements.o \
./example/density/example_density_simulate_default.o \
./example/density/example_density_simulate_gaussian_sphere.o \
./example/density/example_density_simulators.o 

CPP_DEPS += \
./example/density/example_density.d \
./example/density/example_density_connectivity.d \
./example/density/example_density_map.d \
./example/density/example_density_map_cylindrical.d \
./example/density/example_density_mask_3d.d \
./example/density/example_density_protein_agreement_ccc.d \
./example/density/example_density_protein_agreement_likelihood.d \
./example/density/example_density_protein_agreements.d \
./example/density/example_density_simulate_default.d \
./example/density/example_density_simulate_gaussian_sphere.d \
./example/density/example_density_simulators.d 


# Each subdirectory must supply rules for building sources it contributes
example/density/%.o: ../example/density/%.cpp example/density/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



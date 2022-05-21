################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/opti/bcl_opti.cpp \
../source/opti/bcl_opti_ensemble_filter.cpp \
../source/opti/bcl_opti_improvement_type.cpp \
../source/opti/bcl_opti_phase.cpp \
../source/opti/bcl_opti_step_status.cpp \
../source/opti/bcl_opti_template_instantiations.cpp \
../source/opti/bcl_opti_tracker_base.cpp 

OBJS += \
./source/opti/bcl_opti.o \
./source/opti/bcl_opti_ensemble_filter.o \
./source/opti/bcl_opti_improvement_type.o \
./source/opti/bcl_opti_phase.o \
./source/opti/bcl_opti_step_status.o \
./source/opti/bcl_opti_template_instantiations.o \
./source/opti/bcl_opti_tracker_base.o 

CPP_DEPS += \
./source/opti/bcl_opti.d \
./source/opti/bcl_opti_ensemble_filter.d \
./source/opti/bcl_opti_improvement_type.d \
./source/opti/bcl_opti_phase.d \
./source/opti/bcl_opti_step_status.d \
./source/opti/bcl_opti_template_instantiations.d \
./source/opti/bcl_opti_tracker_base.d 


# Each subdirectory must supply rules for building sources it contributes
source/opti/%.o: ../source/opti/%.cpp source/opti/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



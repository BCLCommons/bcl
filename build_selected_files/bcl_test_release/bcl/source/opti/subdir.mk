################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/opti/bcl_opti.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_ensemble_filter.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_improvement_type.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_phase.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_step_status.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_template_instantiations.cpp \
../bcl_test_release/bcl/source/opti/bcl_opti_tracker_base.cpp 

OBJS += \
./bcl_test_release/bcl/source/opti/bcl_opti.o \
./bcl_test_release/bcl/source/opti/bcl_opti_ensemble_filter.o \
./bcl_test_release/bcl/source/opti/bcl_opti_improvement_type.o \
./bcl_test_release/bcl/source/opti/bcl_opti_phase.o \
./bcl_test_release/bcl/source/opti/bcl_opti_step_status.o \
./bcl_test_release/bcl/source/opti/bcl_opti_template_instantiations.o \
./bcl_test_release/bcl/source/opti/bcl_opti_tracker_base.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/opti/bcl_opti.d \
./bcl_test_release/bcl/source/opti/bcl_opti_ensemble_filter.d \
./bcl_test_release/bcl/source/opti/bcl_opti_improvement_type.d \
./bcl_test_release/bcl/source/opti/bcl_opti_phase.d \
./bcl_test_release/bcl/source/opti/bcl_opti_step_status.d \
./bcl_test_release/bcl/source/opti/bcl_opti_template_instantiations.d \
./bcl_test_release/bcl/source/opti/bcl_opti_tracker_base.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/opti/%.o: ../bcl_test_release/bcl/source/opti/%.cpp bcl_test_release/bcl/source/opti/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



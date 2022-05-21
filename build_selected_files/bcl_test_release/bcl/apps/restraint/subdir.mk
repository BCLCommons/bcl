################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/apps/restraint/bcl_app_analyze_spin_label_parameters.cpp \
../bcl_test_release/bcl/apps/restraint/bcl_app_restraint_pofr.cpp \
../bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs.cpp \
../bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs_prep.cpp 

OBJS += \
./bcl_test_release/bcl/apps/restraint/bcl_app_analyze_spin_label_parameters.o \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_pofr.o \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs.o \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs_prep.o 

CPP_DEPS += \
./bcl_test_release/bcl/apps/restraint/bcl_app_analyze_spin_label_parameters.d \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_pofr.d \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs.d \
./bcl_test_release/bcl/apps/restraint/bcl_app_restraint_saxs_prep.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/apps/restraint/%.o: ../bcl_test_release/bcl/apps/restraint/%.cpp bcl_test_release/bcl/apps/restraint/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/restraint/bcl_app_analyze_spin_label_parameters.cpp \
../apps/restraint/bcl_app_restraint_pofr.cpp \
../apps/restraint/bcl_app_restraint_saxs.cpp \
../apps/restraint/bcl_app_restraint_saxs_prep.cpp 

OBJS += \
./apps/restraint/bcl_app_analyze_spin_label_parameters.o \
./apps/restraint/bcl_app_restraint_pofr.o \
./apps/restraint/bcl_app_restraint_saxs.o \
./apps/restraint/bcl_app_restraint_saxs_prep.o 

CPP_DEPS += \
./apps/restraint/bcl_app_analyze_spin_label_parameters.d \
./apps/restraint/bcl_app_restraint_pofr.d \
./apps/restraint/bcl_app_restraint_saxs.d \
./apps/restraint/bcl_app_restraint_saxs_prep.d 


# Each subdirectory must supply rules for building sources it contributes
apps/restraint/%.o: ../apps/restraint/%.cpp apps/restraint/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



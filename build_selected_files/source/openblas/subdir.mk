################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/openblas/bcl_openblas.cpp \
../source/openblas/bcl_openblas_operations.cpp 

OBJS += \
./source/openblas/bcl_openblas.o \
./source/openblas/bcl_openblas_operations.o 

CPP_DEPS += \
./source/openblas/bcl_openblas.d \
./source/openblas/bcl_openblas_operations.d 


# Each subdirectory must supply rules for building sources it contributes
source/openblas/%.o: ../source/openblas/%.cpp source/openblas/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



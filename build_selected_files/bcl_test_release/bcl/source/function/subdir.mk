################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/function/bcl_function.cpp 

OBJS += \
./bcl_test_release/bcl/source/function/bcl_function.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/function/bcl_function.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/function/%.o: ../bcl_test_release/bcl/source/function/%.cpp bcl_test_release/bcl/source/function/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



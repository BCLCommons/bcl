################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/bcl.cpp \
../bcl_test_release/bcl/source/bcl_license.cpp \
../bcl_test_release/bcl/source/bcl_version.cpp 

OBJS += \
./bcl_test_release/bcl/source/bcl.o \
./bcl_test_release/bcl/source/bcl_license.o \
./bcl_test_release/bcl/source/bcl_version.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/bcl.d \
./bcl_test_release/bcl/source/bcl_license.d \
./bcl_test_release/bcl/source/bcl_version.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/%.o: ../bcl_test_release/bcl/source/%.cpp bcl_test_release/bcl/source/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



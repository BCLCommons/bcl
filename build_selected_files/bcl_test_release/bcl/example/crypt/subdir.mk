################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/crypt/example_crypt_blowfish.cpp \
../bcl_test_release/bcl/example/crypt/example_crypt_sha1.cpp 

OBJS += \
./bcl_test_release/bcl/example/crypt/example_crypt_blowfish.o \
./bcl_test_release/bcl/example/crypt/example_crypt_sha1.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/crypt/example_crypt_blowfish.d \
./bcl_test_release/bcl/example/crypt/example_crypt_sha1.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/crypt/%.o: ../bcl_test_release/bcl/example/crypt/%.cpp bcl_test_release/bcl/example/crypt/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



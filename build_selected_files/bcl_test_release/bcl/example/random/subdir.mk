################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/random/example_random.cpp \
../bcl_test_release/bcl/example/random/example_random_histogram_1d_distribution.cpp \
../bcl_test_release/bcl/example/random/example_random_histogram_2d_distribution.cpp \
../bcl_test_release/bcl/example/random/example_random_uniform_distribution.cpp 

OBJS += \
./bcl_test_release/bcl/example/random/example_random.o \
./bcl_test_release/bcl/example/random/example_random_histogram_1d_distribution.o \
./bcl_test_release/bcl/example/random/example_random_histogram_2d_distribution.o \
./bcl_test_release/bcl/example/random/example_random_uniform_distribution.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/random/example_random.d \
./bcl_test_release/bcl/example/random/example_random_histogram_1d_distribution.d \
./bcl_test_release/bcl/example/random/example_random_histogram_2d_distribution.d \
./bcl_test_release/bcl/example/random/example_random_uniform_distribution.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/random/%.o: ../bcl_test_release/bcl/example/random/%.cpp bcl_test_release/bcl/example/random/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



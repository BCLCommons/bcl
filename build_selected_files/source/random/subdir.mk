################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/random/bcl_random.cpp \
../source/random/bcl_random_distribution_interface.cpp \
../source/random/bcl_random_histogram_1d_distribution.cpp \
../source/random/bcl_random_histogram_2d_distribution.cpp \
../source/random/bcl_random_uniform_distribution.cpp 

OBJS += \
./source/random/bcl_random.o \
./source/random/bcl_random_distribution_interface.o \
./source/random/bcl_random_histogram_1d_distribution.o \
./source/random/bcl_random_histogram_2d_distribution.o \
./source/random/bcl_random_uniform_distribution.o 

CPP_DEPS += \
./source/random/bcl_random.d \
./source/random/bcl_random_distribution_interface.d \
./source/random/bcl_random_histogram_1d_distribution.d \
./source/random/bcl_random_histogram_2d_distribution.d \
./source/random/bcl_random_uniform_distribution.d 


# Each subdirectory must supply rules for building sources it contributes
source/random/%.o: ../source/random/%.cpp source/random/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/bzip2/bcl_bzip2.cpp \
../source/bzip2/bcl_bzip2_stream_buffer.cpp 

OBJS += \
./source/bzip2/bcl_bzip2.o \
./source/bzip2/bcl_bzip2_stream_buffer.o 

CPP_DEPS += \
./source/bzip2/bcl_bzip2.d \
./source/bzip2/bcl_bzip2_stream_buffer.d 


# Each subdirectory must supply rules for building sources it contributes
source/bzip2/%.o: ../source/bzip2/%.cpp source/bzip2/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



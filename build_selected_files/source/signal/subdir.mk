################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/signal/bcl_signal.cpp \
../source/signal/bcl_signal_slots.cpp 

OBJS += \
./source/signal/bcl_signal.o \
./source/signal/bcl_signal_slots.o 

CPP_DEPS += \
./source/signal/bcl_signal.d \
./source/signal/bcl_signal_slots.d 


# Each subdirectory must supply rules for building sources it contributes
source/signal/%.o: ../source/signal/%.cpp source/signal/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



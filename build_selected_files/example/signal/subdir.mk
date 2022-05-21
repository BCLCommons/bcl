################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/signal/example_signal.cpp \
../example/signal/example_signal_connection.cpp \
../example/signal/example_signal_signal.cpp \
../example/signal/example_signal_slots.cpp 

OBJS += \
./example/signal/example_signal.o \
./example/signal/example_signal_connection.o \
./example/signal/example_signal_signal.o \
./example/signal/example_signal_slots.o 

CPP_DEPS += \
./example/signal/example_signal.d \
./example/signal/example_signal_connection.d \
./example/signal/example_signal_signal.d \
./example/signal/example_signal_slots.d 


# Each subdirectory must supply rules for building sources it contributes
example/signal/%.o: ../example/signal/%.cpp example/signal/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



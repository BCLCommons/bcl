################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/function/example_function.cpp \
../example/function/example_function_binary_cached.cpp \
../example/function/example_function_member.cpp \
../example/function/example_function_member_const.cpp \
../example/function/example_function_member_unary.cpp \
../example/function/example_function_member_unary_const.cpp \
../example/function/example_function_unary_cached.cpp 

OBJS += \
./example/function/example_function.o \
./example/function/example_function_binary_cached.o \
./example/function/example_function_member.o \
./example/function/example_function_member_const.o \
./example/function/example_function_member_unary.o \
./example/function/example_function_member_unary_const.o \
./example/function/example_function_unary_cached.o 

CPP_DEPS += \
./example/function/example_function.d \
./example/function/example_function_binary_cached.d \
./example/function/example_function_member.d \
./example/function/example_function_member_const.d \
./example/function/example_function_member_unary.d \
./example/function/example_function_member_unary_const.d \
./example/function/example_function_unary_cached.d 


# Each subdirectory must supply rules for building sources it contributes
example/function/%.o: ../example/function/%.cpp example/function/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



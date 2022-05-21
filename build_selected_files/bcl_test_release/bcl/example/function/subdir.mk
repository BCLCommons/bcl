################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/function/example_function.cpp \
../bcl_test_release/bcl/example/function/example_function_binary_cached.cpp \
../bcl_test_release/bcl/example/function/example_function_member.cpp \
../bcl_test_release/bcl/example/function/example_function_member_const.cpp \
../bcl_test_release/bcl/example/function/example_function_member_unary.cpp \
../bcl_test_release/bcl/example/function/example_function_member_unary_const.cpp \
../bcl_test_release/bcl/example/function/example_function_unary_cached.cpp 

OBJS += \
./bcl_test_release/bcl/example/function/example_function.o \
./bcl_test_release/bcl/example/function/example_function_binary_cached.o \
./bcl_test_release/bcl/example/function/example_function_member.o \
./bcl_test_release/bcl/example/function/example_function_member_const.o \
./bcl_test_release/bcl/example/function/example_function_member_unary.o \
./bcl_test_release/bcl/example/function/example_function_member_unary_const.o \
./bcl_test_release/bcl/example/function/example_function_unary_cached.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/function/example_function.d \
./bcl_test_release/bcl/example/function/example_function_binary_cached.d \
./bcl_test_release/bcl/example/function/example_function_member.d \
./bcl_test_release/bcl/example/function/example_function_member_const.d \
./bcl_test_release/bcl/example/function/example_function_member_unary.d \
./bcl_test_release/bcl/example/function/example_function_member_unary_const.d \
./bcl_test_release/bcl/example/function/example_function_unary_cached.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/function/%.o: ../bcl_test_release/bcl/example/function/%.cpp bcl_test_release/bcl/example/function/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



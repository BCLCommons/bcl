################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/iterate/example_iterate.cpp \
../bcl_test_release/bcl/example/iterate/example_iterate_generic.cpp \
../bcl_test_release/bcl/example/iterate/example_iterate_reflecting.cpp \
../bcl_test_release/bcl/example/iterate/example_iterate_with_range.cpp 

OBJS += \
./bcl_test_release/bcl/example/iterate/example_iterate.o \
./bcl_test_release/bcl/example/iterate/example_iterate_generic.o \
./bcl_test_release/bcl/example/iterate/example_iterate_reflecting.o \
./bcl_test_release/bcl/example/iterate/example_iterate_with_range.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/iterate/example_iterate.d \
./bcl_test_release/bcl/example/iterate/example_iterate_generic.d \
./bcl_test_release/bcl/example/iterate/example_iterate_reflecting.d \
./bcl_test_release/bcl/example/iterate/example_iterate_with_range.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/iterate/%.o: ../bcl_test_release/bcl/example/iterate/%.cpp bcl_test_release/bcl/example/iterate/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



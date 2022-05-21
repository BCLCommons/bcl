################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/type/example_type.cpp \
../bcl_test_release/bcl/example/type/example_type_chooser.cpp \
../bcl_test_release/bcl/example/type/example_type_compare.cpp \
../bcl_test_release/bcl/example/type/example_type_enable_if.cpp \
../bcl_test_release/bcl/example/type/example_type_is_a.cpp \
../bcl_test_release/bcl/example/type/example_type_is_derived_from.cpp \
../bcl_test_release/bcl/example/type/example_type_is_map.cpp \
../bcl_test_release/bcl/example/type/example_type_is_sequence.cpp \
../bcl_test_release/bcl/example/type/example_type_remove_const_ref.cpp 

OBJS += \
./bcl_test_release/bcl/example/type/example_type.o \
./bcl_test_release/bcl/example/type/example_type_chooser.o \
./bcl_test_release/bcl/example/type/example_type_compare.o \
./bcl_test_release/bcl/example/type/example_type_enable_if.o \
./bcl_test_release/bcl/example/type/example_type_is_a.o \
./bcl_test_release/bcl/example/type/example_type_is_derived_from.o \
./bcl_test_release/bcl/example/type/example_type_is_map.o \
./bcl_test_release/bcl/example/type/example_type_is_sequence.o \
./bcl_test_release/bcl/example/type/example_type_remove_const_ref.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/type/example_type.d \
./bcl_test_release/bcl/example/type/example_type_chooser.d \
./bcl_test_release/bcl/example/type/example_type_compare.d \
./bcl_test_release/bcl/example/type/example_type_enable_if.d \
./bcl_test_release/bcl/example/type/example_type_is_a.d \
./bcl_test_release/bcl/example/type/example_type_is_derived_from.d \
./bcl_test_release/bcl/example/type/example_type_is_map.d \
./bcl_test_release/bcl/example/type/example_type_is_sequence.d \
./bcl_test_release/bcl/example/type/example_type_remove_const_ref.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/type/%.o: ../bcl_test_release/bcl/example/type/%.cpp bcl_test_release/bcl/example/type/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



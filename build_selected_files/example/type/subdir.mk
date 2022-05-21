################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/type/example_type.cpp \
../example/type/example_type_chooser.cpp \
../example/type/example_type_compare.cpp \
../example/type/example_type_enable_if.cpp \
../example/type/example_type_is_a.cpp \
../example/type/example_type_is_derived_from.cpp \
../example/type/example_type_is_map.cpp \
../example/type/example_type_is_sequence.cpp \
../example/type/example_type_remove_const_ref.cpp 

OBJS += \
./example/type/example_type.o \
./example/type/example_type_chooser.o \
./example/type/example_type_compare.o \
./example/type/example_type_enable_if.o \
./example/type/example_type_is_a.o \
./example/type/example_type_is_derived_from.o \
./example/type/example_type_is_map.o \
./example/type/example_type_is_sequence.o \
./example/type/example_type_remove_const_ref.o 

CPP_DEPS += \
./example/type/example_type.d \
./example/type/example_type_chooser.d \
./example/type/example_type_compare.d \
./example/type/example_type_enable_if.d \
./example/type/example_type_is_a.d \
./example/type/example_type_is_derived_from.d \
./example/type/example_type_is_map.d \
./example/type/example_type_is_sequence.d \
./example/type/example_type_remove_const_ref.d 


# Each subdirectory must supply rules for building sources it contributes
example/type/%.o: ../example/type/%.cpp example/type/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



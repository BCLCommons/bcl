################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/storage/bcl_storage.cpp \
../bcl_test_release/bcl/source/storage/bcl_storage_table.cpp \
../bcl_test_release/bcl/source/storage/bcl_storage_table_header.cpp \
../bcl_test_release/bcl/source/storage/bcl_storage_template_instantiations.cpp \
../bcl_test_release/bcl/source/storage/bcl_storage_vector.cpp 

OBJS += \
./bcl_test_release/bcl/source/storage/bcl_storage.o \
./bcl_test_release/bcl/source/storage/bcl_storage_table.o \
./bcl_test_release/bcl/source/storage/bcl_storage_table_header.o \
./bcl_test_release/bcl/source/storage/bcl_storage_template_instantiations.o \
./bcl_test_release/bcl/source/storage/bcl_storage_vector.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/storage/bcl_storage.d \
./bcl_test_release/bcl/source/storage/bcl_storage_table.d \
./bcl_test_release/bcl/source/storage/bcl_storage_table_header.d \
./bcl_test_release/bcl/source/storage/bcl_storage_template_instantiations.d \
./bcl_test_release/bcl/source/storage/bcl_storage_vector.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/storage/%.o: ../bcl_test_release/bcl/source/storage/%.cpp bcl_test_release/bcl/source/storage/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



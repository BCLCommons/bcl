################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/storage/bcl_storage.cpp \
../source/storage/bcl_storage_table.cpp \
../source/storage/bcl_storage_table_header.cpp \
../source/storage/bcl_storage_template_instantiations.cpp \
../source/storage/bcl_storage_vector.cpp 

OBJS += \
./source/storage/bcl_storage.o \
./source/storage/bcl_storage_table.o \
./source/storage/bcl_storage_table_header.o \
./source/storage/bcl_storage_template_instantiations.o \
./source/storage/bcl_storage_vector.o 

CPP_DEPS += \
./source/storage/bcl_storage.d \
./source/storage/bcl_storage_table.d \
./source/storage/bcl_storage_table_header.d \
./source/storage/bcl_storage_template_instantiations.d \
./source/storage/bcl_storage_vector.d 


# Each subdirectory must supply rules for building sources it contributes
source/storage/%.o: ../source/storage/%.cpp source/storage/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



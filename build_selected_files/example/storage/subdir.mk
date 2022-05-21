################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/storage/example_storage.cpp \
../example/storage/example_storage_hash_map.cpp \
../example/storage/example_storage_list.cpp \
../example/storage/example_storage_map.cpp \
../example/storage/example_storage_object_nd_hash_map.cpp \
../example/storage/example_storage_pair.cpp \
../example/storage/example_storage_row.cpp \
../example/storage/example_storage_set.cpp \
../example/storage/example_storage_symmetric_matrix.cpp \
../example/storage/example_storage_table.cpp \
../example/storage/example_storage_table_header.cpp \
../example/storage/example_storage_triplet.cpp \
../example/storage/example_storage_vector.cpp \
../example/storage/example_storage_vector_nd.cpp 

OBJS += \
./example/storage/example_storage.o \
./example/storage/example_storage_hash_map.o \
./example/storage/example_storage_list.o \
./example/storage/example_storage_map.o \
./example/storage/example_storage_object_nd_hash_map.o \
./example/storage/example_storage_pair.o \
./example/storage/example_storage_row.o \
./example/storage/example_storage_set.o \
./example/storage/example_storage_symmetric_matrix.o \
./example/storage/example_storage_table.o \
./example/storage/example_storage_table_header.o \
./example/storage/example_storage_triplet.o \
./example/storage/example_storage_vector.o \
./example/storage/example_storage_vector_nd.o 

CPP_DEPS += \
./example/storage/example_storage.d \
./example/storage/example_storage_hash_map.d \
./example/storage/example_storage_list.d \
./example/storage/example_storage_map.d \
./example/storage/example_storage_object_nd_hash_map.d \
./example/storage/example_storage_pair.d \
./example/storage/example_storage_row.d \
./example/storage/example_storage_set.d \
./example/storage/example_storage_symmetric_matrix.d \
./example/storage/example_storage_table.d \
./example/storage/example_storage_table_header.d \
./example/storage/example_storage_triplet.d \
./example/storage/example_storage_vector.d \
./example/storage/example_storage_vector_nd.d 


# Each subdirectory must supply rules for building sources it contributes
example/storage/%.o: ../example/storage/%.cpp example/storage/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



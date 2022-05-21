################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/storage/example_storage.cpp \
../bcl_test_release/bcl/example/storage/example_storage_hash_map.cpp \
../bcl_test_release/bcl/example/storage/example_storage_list.cpp \
../bcl_test_release/bcl/example/storage/example_storage_map.cpp \
../bcl_test_release/bcl/example/storage/example_storage_object_nd_hash_map.cpp \
../bcl_test_release/bcl/example/storage/example_storage_pair.cpp \
../bcl_test_release/bcl/example/storage/example_storage_row.cpp \
../bcl_test_release/bcl/example/storage/example_storage_set.cpp \
../bcl_test_release/bcl/example/storage/example_storage_symmetric_matrix.cpp \
../bcl_test_release/bcl/example/storage/example_storage_table.cpp \
../bcl_test_release/bcl/example/storage/example_storage_table_header.cpp \
../bcl_test_release/bcl/example/storage/example_storage_triplet.cpp \
../bcl_test_release/bcl/example/storage/example_storage_vector.cpp \
../bcl_test_release/bcl/example/storage/example_storage_vector_nd.cpp 

OBJS += \
./bcl_test_release/bcl/example/storage/example_storage.o \
./bcl_test_release/bcl/example/storage/example_storage_hash_map.o \
./bcl_test_release/bcl/example/storage/example_storage_list.o \
./bcl_test_release/bcl/example/storage/example_storage_map.o \
./bcl_test_release/bcl/example/storage/example_storage_object_nd_hash_map.o \
./bcl_test_release/bcl/example/storage/example_storage_pair.o \
./bcl_test_release/bcl/example/storage/example_storage_row.o \
./bcl_test_release/bcl/example/storage/example_storage_set.o \
./bcl_test_release/bcl/example/storage/example_storage_symmetric_matrix.o \
./bcl_test_release/bcl/example/storage/example_storage_table.o \
./bcl_test_release/bcl/example/storage/example_storage_table_header.o \
./bcl_test_release/bcl/example/storage/example_storage_triplet.o \
./bcl_test_release/bcl/example/storage/example_storage_vector.o \
./bcl_test_release/bcl/example/storage/example_storage_vector_nd.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/storage/example_storage.d \
./bcl_test_release/bcl/example/storage/example_storage_hash_map.d \
./bcl_test_release/bcl/example/storage/example_storage_list.d \
./bcl_test_release/bcl/example/storage/example_storage_map.d \
./bcl_test_release/bcl/example/storage/example_storage_object_nd_hash_map.d \
./bcl_test_release/bcl/example/storage/example_storage_pair.d \
./bcl_test_release/bcl/example/storage/example_storage_row.d \
./bcl_test_release/bcl/example/storage/example_storage_set.d \
./bcl_test_release/bcl/example/storage/example_storage_symmetric_matrix.d \
./bcl_test_release/bcl/example/storage/example_storage_table.d \
./bcl_test_release/bcl/example/storage/example_storage_table_header.d \
./bcl_test_release/bcl/example/storage/example_storage_triplet.d \
./bcl_test_release/bcl/example/storage/example_storage_vector.d \
./bcl_test_release/bcl/example/storage/example_storage_vector_nd.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/storage/%.o: ../bcl_test_release/bcl/example/storage/%.cpp bcl_test_release/bcl/example/storage/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



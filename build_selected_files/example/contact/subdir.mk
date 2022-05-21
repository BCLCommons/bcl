################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/contact/example_contact.cpp \
../example/contact/example_contact_aa_correlation_from_file.cpp \
../example/contact/example_contact_calculate_correlations_mi.cpp \
../example/contact/example_contact_calculate_correlations_sm.cpp \
../example/contact/example_contact_correlation_matrix.cpp \
../example/contact/example_contact_correlation_storage_file.cpp \
../example/contact/example_contact_data.cpp \
../example/contact/example_contact_map.cpp \
../example/contact/example_contact_order.cpp \
../example/contact/example_contact_prediction_map.cpp \
../example/contact/example_contact_recovery.cpp \
../example/contact/example_contact_sse_prediction_map.cpp \
../example/contact/example_contact_statistics.cpp \
../example/contact/example_contact_type_data.cpp \
../example/contact/example_contact_types.cpp 

OBJS += \
./example/contact/example_contact.o \
./example/contact/example_contact_aa_correlation_from_file.o \
./example/contact/example_contact_calculate_correlations_mi.o \
./example/contact/example_contact_calculate_correlations_sm.o \
./example/contact/example_contact_correlation_matrix.o \
./example/contact/example_contact_correlation_storage_file.o \
./example/contact/example_contact_data.o \
./example/contact/example_contact_map.o \
./example/contact/example_contact_order.o \
./example/contact/example_contact_prediction_map.o \
./example/contact/example_contact_recovery.o \
./example/contact/example_contact_sse_prediction_map.o \
./example/contact/example_contact_statistics.o \
./example/contact/example_contact_type_data.o \
./example/contact/example_contact_types.o 

CPP_DEPS += \
./example/contact/example_contact.d \
./example/contact/example_contact_aa_correlation_from_file.d \
./example/contact/example_contact_calculate_correlations_mi.d \
./example/contact/example_contact_calculate_correlations_sm.d \
./example/contact/example_contact_correlation_matrix.d \
./example/contact/example_contact_correlation_storage_file.d \
./example/contact/example_contact_data.d \
./example/contact/example_contact_map.d \
./example/contact/example_contact_order.d \
./example/contact/example_contact_prediction_map.d \
./example/contact/example_contact_recovery.d \
./example/contact/example_contact_sse_prediction_map.d \
./example/contact/example_contact_statistics.d \
./example/contact/example_contact_type_data.d \
./example/contact/example_contact_types.d 


# Each subdirectory must supply rules for building sources it contributes
example/contact/%.o: ../example/contact/%.cpp example/contact/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



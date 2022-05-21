################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/contact/bcl_contact.cpp \
../source/contact/bcl_contact_aa_correlation_from_file.cpp \
../source/contact/bcl_contact_ann.cpp \
../source/contact/bcl_contact_calculate_correlations_mi.cpp \
../source/contact/bcl_contact_calculate_correlations_sm.cpp \
../source/contact/bcl_contact_correlation_matrix.cpp \
../source/contact/bcl_contact_correlation_storage_file.cpp \
../source/contact/bcl_contact_data.cpp \
../source/contact/bcl_contact_map.cpp \
../source/contact/bcl_contact_order.cpp \
../source/contact/bcl_contact_prediction_map.cpp \
../source/contact/bcl_contact_recovery.cpp \
../source/contact/bcl_contact_sse_prediction_map.cpp \
../source/contact/bcl_contact_statistics.cpp \
../source/contact/bcl_contact_type_data.cpp \
../source/contact/bcl_contact_types.cpp 

OBJS += \
./source/contact/bcl_contact.o \
./source/contact/bcl_contact_aa_correlation_from_file.o \
./source/contact/bcl_contact_ann.o \
./source/contact/bcl_contact_calculate_correlations_mi.o \
./source/contact/bcl_contact_calculate_correlations_sm.o \
./source/contact/bcl_contact_correlation_matrix.o \
./source/contact/bcl_contact_correlation_storage_file.o \
./source/contact/bcl_contact_data.o \
./source/contact/bcl_contact_map.o \
./source/contact/bcl_contact_order.o \
./source/contact/bcl_contact_prediction_map.o \
./source/contact/bcl_contact_recovery.o \
./source/contact/bcl_contact_sse_prediction_map.o \
./source/contact/bcl_contact_statistics.o \
./source/contact/bcl_contact_type_data.o \
./source/contact/bcl_contact_types.o 

CPP_DEPS += \
./source/contact/bcl_contact.d \
./source/contact/bcl_contact_aa_correlation_from_file.d \
./source/contact/bcl_contact_ann.d \
./source/contact/bcl_contact_calculate_correlations_mi.d \
./source/contact/bcl_contact_calculate_correlations_sm.d \
./source/contact/bcl_contact_correlation_matrix.d \
./source/contact/bcl_contact_correlation_storage_file.d \
./source/contact/bcl_contact_data.d \
./source/contact/bcl_contact_map.d \
./source/contact/bcl_contact_order.d \
./source/contact/bcl_contact_prediction_map.d \
./source/contact/bcl_contact_recovery.d \
./source/contact/bcl_contact_sse_prediction_map.d \
./source/contact/bcl_contact_statistics.d \
./source/contact/bcl_contact_type_data.d \
./source/contact/bcl_contact_types.d 


# Each subdirectory must supply rules for building sources it contributes
source/contact/%.o: ../source/contact/%.cpp source/contact/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



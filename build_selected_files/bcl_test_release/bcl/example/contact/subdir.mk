################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/contact/example_contact.cpp \
../bcl_test_release/bcl/example/contact/example_contact_aa_correlation_from_file.cpp \
../bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_mi.cpp \
../bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_sm.cpp \
../bcl_test_release/bcl/example/contact/example_contact_correlation_matrix.cpp \
../bcl_test_release/bcl/example/contact/example_contact_correlation_storage_file.cpp \
../bcl_test_release/bcl/example/contact/example_contact_data.cpp \
../bcl_test_release/bcl/example/contact/example_contact_map.cpp \
../bcl_test_release/bcl/example/contact/example_contact_order.cpp \
../bcl_test_release/bcl/example/contact/example_contact_prediction_map.cpp \
../bcl_test_release/bcl/example/contact/example_contact_recovery.cpp \
../bcl_test_release/bcl/example/contact/example_contact_sse_prediction_map.cpp \
../bcl_test_release/bcl/example/contact/example_contact_statistics.cpp \
../bcl_test_release/bcl/example/contact/example_contact_type_data.cpp \
../bcl_test_release/bcl/example/contact/example_contact_types.cpp 

OBJS += \
./bcl_test_release/bcl/example/contact/example_contact.o \
./bcl_test_release/bcl/example/contact/example_contact_aa_correlation_from_file.o \
./bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_mi.o \
./bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_sm.o \
./bcl_test_release/bcl/example/contact/example_contact_correlation_matrix.o \
./bcl_test_release/bcl/example/contact/example_contact_correlation_storage_file.o \
./bcl_test_release/bcl/example/contact/example_contact_data.o \
./bcl_test_release/bcl/example/contact/example_contact_map.o \
./bcl_test_release/bcl/example/contact/example_contact_order.o \
./bcl_test_release/bcl/example/contact/example_contact_prediction_map.o \
./bcl_test_release/bcl/example/contact/example_contact_recovery.o \
./bcl_test_release/bcl/example/contact/example_contact_sse_prediction_map.o \
./bcl_test_release/bcl/example/contact/example_contact_statistics.o \
./bcl_test_release/bcl/example/contact/example_contact_type_data.o \
./bcl_test_release/bcl/example/contact/example_contact_types.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/contact/example_contact.d \
./bcl_test_release/bcl/example/contact/example_contact_aa_correlation_from_file.d \
./bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_mi.d \
./bcl_test_release/bcl/example/contact/example_contact_calculate_correlations_sm.d \
./bcl_test_release/bcl/example/contact/example_contact_correlation_matrix.d \
./bcl_test_release/bcl/example/contact/example_contact_correlation_storage_file.d \
./bcl_test_release/bcl/example/contact/example_contact_data.d \
./bcl_test_release/bcl/example/contact/example_contact_map.d \
./bcl_test_release/bcl/example/contact/example_contact_order.d \
./bcl_test_release/bcl/example/contact/example_contact_prediction_map.d \
./bcl_test_release/bcl/example/contact/example_contact_recovery.d \
./bcl_test_release/bcl/example/contact/example_contact_sse_prediction_map.d \
./bcl_test_release/bcl/example/contact/example_contact_statistics.d \
./bcl_test_release/bcl/example/contact/example_contact_type_data.d \
./bcl_test_release/bcl/example/contact/example_contact_types.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/contact/%.o: ../bcl_test_release/bcl/example/contact/%.cpp bcl_test_release/bcl/example/contact/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



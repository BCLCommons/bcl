################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/apps/internal/util/bcl_app_cpu_benchmark.cpp \
../bcl_test_release/bcl/apps/internal/util/bcl_app_generate_license_file.cpp \
../bcl_test_release/bcl/apps/internal/util/bcl_app_minimize_score_weight_set.cpp \
../bcl_test_release/bcl/apps/internal/util/bcl_app_write_app_web_text.cpp \
../bcl_test_release/bcl/apps/internal/util/bcl_app_write_bcl_object.cpp 

OBJS += \
./bcl_test_release/bcl/apps/internal/util/bcl_app_cpu_benchmark.o \
./bcl_test_release/bcl/apps/internal/util/bcl_app_generate_license_file.o \
./bcl_test_release/bcl/apps/internal/util/bcl_app_minimize_score_weight_set.o \
./bcl_test_release/bcl/apps/internal/util/bcl_app_write_app_web_text.o \
./bcl_test_release/bcl/apps/internal/util/bcl_app_write_bcl_object.o 

CPP_DEPS += \
./bcl_test_release/bcl/apps/internal/util/bcl_app_cpu_benchmark.d \
./bcl_test_release/bcl/apps/internal/util/bcl_app_generate_license_file.d \
./bcl_test_release/bcl/apps/internal/util/bcl_app_minimize_score_weight_set.d \
./bcl_test_release/bcl/apps/internal/util/bcl_app_write_app_web_text.d \
./bcl_test_release/bcl/apps/internal/util/bcl_app_write_bcl_object.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/apps/internal/util/%.o: ../bcl_test_release/bcl/apps/internal/util/%.cpp bcl_test_release/bcl/apps/internal/util/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



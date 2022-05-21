################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/apps/model/example_app_model_compute_statistics.cpp \
../bcl_test_release/bcl/example/apps/model/example_app_model_prediction_merge.cpp \
../bcl_test_release/bcl/example/apps/model/example_app_model_test.cpp \
../bcl_test_release/bcl/example/apps/model/example_app_model_train.cpp 

OBJS += \
./bcl_test_release/bcl/example/apps/model/example_app_model_compute_statistics.o \
./bcl_test_release/bcl/example/apps/model/example_app_model_prediction_merge.o \
./bcl_test_release/bcl/example/apps/model/example_app_model_test.o \
./bcl_test_release/bcl/example/apps/model/example_app_model_train.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/apps/model/example_app_model_compute_statistics.d \
./bcl_test_release/bcl/example/apps/model/example_app_model_prediction_merge.d \
./bcl_test_release/bcl/example/apps/model/example_app_model_test.d \
./bcl_test_release/bcl/example/apps/model/example_app_model_train.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/apps/model/%.o: ../bcl_test_release/bcl/example/apps/model/%.cpp bcl_test_release/bcl/example/apps/model/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



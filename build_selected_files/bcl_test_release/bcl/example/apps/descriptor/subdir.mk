################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_analyze.cpp \
../bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_generate_dataset.cpp \
../bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_refine_by_score.cpp \
../bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_score_dataset.cpp \
../bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_sequential_feature_selection.cpp 

OBJS += \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_analyze.o \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_generate_dataset.o \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_refine_by_score.o \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_score_dataset.o \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_sequential_feature_selection.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_analyze.d \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_generate_dataset.d \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_refine_by_score.d \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_score_dataset.d \
./bcl_test_release/bcl/example/apps/descriptor/example_app_descriptor_sequential_feature_selection.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/apps/descriptor/%.o: ../bcl_test_release/bcl/example/apps/descriptor/%.cpp bcl_test_release/bcl/example/apps/descriptor/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



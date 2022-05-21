################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/bcl_app_examples.cpp \
../bcl_test_release/bcl/example/example.cpp \
../bcl_test_release/bcl/example/example_application_example_helper.cpp \
../bcl_test_release/bcl/example/example_bcl.cpp \
../bcl_test_release/bcl/example/example_bcl_license.cpp \
../bcl_test_release/bcl/example/example_bcl_version.cpp \
../bcl_test_release/bcl/example/example_interface.cpp \
../bcl_test_release/bcl/example/example_proteins.cpp 

OBJS += \
./bcl_test_release/bcl/example/bcl_app_examples.o \
./bcl_test_release/bcl/example/example.o \
./bcl_test_release/bcl/example/example_application_example_helper.o \
./bcl_test_release/bcl/example/example_bcl.o \
./bcl_test_release/bcl/example/example_bcl_license.o \
./bcl_test_release/bcl/example/example_bcl_version.o \
./bcl_test_release/bcl/example/example_interface.o \
./bcl_test_release/bcl/example/example_proteins.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/bcl_app_examples.d \
./bcl_test_release/bcl/example/example.d \
./bcl_test_release/bcl/example/example_application_example_helper.d \
./bcl_test_release/bcl/example/example_bcl.d \
./bcl_test_release/bcl/example/example_bcl_license.d \
./bcl_test_release/bcl/example/example_bcl_version.d \
./bcl_test_release/bcl/example/example_interface.d \
./bcl_test_release/bcl/example/example_proteins.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/%.o: ../bcl_test_release/bcl/example/%.cpp bcl_test_release/bcl/example/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



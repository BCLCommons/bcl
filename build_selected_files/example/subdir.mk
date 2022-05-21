################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/bcl_app_examples.cpp \
../example/example.cpp \
../example/example_application_example_helper.cpp \
../example/example_bcl.cpp \
../example/example_bcl_license.cpp \
../example/example_bcl_version.cpp \
../example/example_interface.cpp \
../example/example_proteins.cpp 

OBJS += \
./example/bcl_app_examples.o \
./example/example.o \
./example/example_application_example_helper.o \
./example/example_bcl.o \
./example/example_bcl_license.o \
./example/example_bcl_version.o \
./example/example_interface.o \
./example/example_proteins.o 

CPP_DEPS += \
./example/bcl_app_examples.d \
./example/example.d \
./example/example_application_example_helper.d \
./example/example_bcl.d \
./example/example_bcl_license.d \
./example/example_bcl_version.d \
./example/example_interface.d \
./example/example_proteins.d 


# Each subdirectory must supply rules for building sources it contributes
example/%.o: ../example/%.cpp example/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



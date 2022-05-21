################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/app/bcl_app.cpp \
../bcl_test_release/bcl/source/app/bcl_app_apps.cpp \
../bcl_test_release/bcl/source/app/bcl_app_group_handler.cpp \
../bcl_test_release/bcl/source/app/bcl_app_groups.cpp \
../bcl_test_release/bcl/source/app/bcl_app_help.cpp \
../bcl_test_release/bcl/source/app/bcl_app_interface.cpp 

OBJS += \
./bcl_test_release/bcl/source/app/bcl_app.o \
./bcl_test_release/bcl/source/app/bcl_app_apps.o \
./bcl_test_release/bcl/source/app/bcl_app_group_handler.o \
./bcl_test_release/bcl/source/app/bcl_app_groups.o \
./bcl_test_release/bcl/source/app/bcl_app_help.o \
./bcl_test_release/bcl/source/app/bcl_app_interface.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/app/bcl_app.d \
./bcl_test_release/bcl/source/app/bcl_app_apps.d \
./bcl_test_release/bcl/source/app/bcl_app_group_handler.d \
./bcl_test_release/bcl/source/app/bcl_app_groups.d \
./bcl_test_release/bcl/source/app/bcl_app_help.d \
./bcl_test_release/bcl/source/app/bcl_app_interface.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/app/%.o: ../bcl_test_release/bcl/source/app/%.cpp bcl_test_release/bcl/source/app/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



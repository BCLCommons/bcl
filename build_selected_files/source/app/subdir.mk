################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/app/bcl_app.cpp \
../source/app/bcl_app_apps.cpp \
../source/app/bcl_app_group_handler.cpp \
../source/app/bcl_app_groups.cpp \
../source/app/bcl_app_help.cpp \
../source/app/bcl_app_interface.cpp 

OBJS += \
./source/app/bcl_app.o \
./source/app/bcl_app_apps.o \
./source/app/bcl_app_group_handler.o \
./source/app/bcl_app_groups.o \
./source/app/bcl_app_help.o \
./source/app/bcl_app_interface.o 

CPP_DEPS += \
./source/app/bcl_app.d \
./source/app/bcl_app_apps.d \
./source/app/bcl_app_group_handler.d \
./source/app/bcl_app_groups.d \
./source/app/bcl_app_help.d \
./source/app/bcl_app_interface.d 


# Each subdirectory must supply rules for building sources it contributes
source/app/%.o: ../source/app/%.cpp source/app/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/command/bcl_command.cpp \
../source/command/bcl_command_app_default_flags.cpp \
../source/command/bcl_command_command.cpp \
../source/command/bcl_command_command_line_writer.cpp \
../source/command/bcl_command_command_state.cpp \
../source/command/bcl_command_default_flag_types.cpp \
../source/command/bcl_command_flag_dynamic.cpp \
../source/command/bcl_command_flag_separator.cpp \
../source/command/bcl_command_flag_static.cpp \
../source/command/bcl_command_flag_static_and_dynamic.cpp \
../source/command/bcl_command_guesser.cpp \
../source/command/bcl_command_parameter.cpp \
../source/command/bcl_command_parameter_check_allowed.cpp \
../source/command/bcl_command_parameter_check_allowed_non_const.cpp \
../source/command/bcl_command_parameter_check_default.cpp \
../source/command/bcl_command_parameter_check_extension.cpp \
../source/command/bcl_command_parameter_check_extensions_file_existence.cpp \
../source/command/bcl_command_parameter_check_file_existence.cpp \
../source/command/bcl_command_parameter_check_file_in_search_path.cpp \
../source/command/bcl_command_parameter_check_interface.cpp \
../source/command/bcl_command_parameter_check_or.cpp \
../source/command/bcl_command_parameter_check_ranged.cpp \
../source/command/bcl_command_parameter_check_serializable.cpp 

OBJS += \
./source/command/bcl_command.o \
./source/command/bcl_command_app_default_flags.o \
./source/command/bcl_command_command.o \
./source/command/bcl_command_command_line_writer.o \
./source/command/bcl_command_command_state.o \
./source/command/bcl_command_default_flag_types.o \
./source/command/bcl_command_flag_dynamic.o \
./source/command/bcl_command_flag_separator.o \
./source/command/bcl_command_flag_static.o \
./source/command/bcl_command_flag_static_and_dynamic.o \
./source/command/bcl_command_guesser.o \
./source/command/bcl_command_parameter.o \
./source/command/bcl_command_parameter_check_allowed.o \
./source/command/bcl_command_parameter_check_allowed_non_const.o \
./source/command/bcl_command_parameter_check_default.o \
./source/command/bcl_command_parameter_check_extension.o \
./source/command/bcl_command_parameter_check_extensions_file_existence.o \
./source/command/bcl_command_parameter_check_file_existence.o \
./source/command/bcl_command_parameter_check_file_in_search_path.o \
./source/command/bcl_command_parameter_check_interface.o \
./source/command/bcl_command_parameter_check_or.o \
./source/command/bcl_command_parameter_check_ranged.o \
./source/command/bcl_command_parameter_check_serializable.o 

CPP_DEPS += \
./source/command/bcl_command.d \
./source/command/bcl_command_app_default_flags.d \
./source/command/bcl_command_command.d \
./source/command/bcl_command_command_line_writer.d \
./source/command/bcl_command_command_state.d \
./source/command/bcl_command_default_flag_types.d \
./source/command/bcl_command_flag_dynamic.d \
./source/command/bcl_command_flag_separator.d \
./source/command/bcl_command_flag_static.d \
./source/command/bcl_command_flag_static_and_dynamic.d \
./source/command/bcl_command_guesser.d \
./source/command/bcl_command_parameter.d \
./source/command/bcl_command_parameter_check_allowed.d \
./source/command/bcl_command_parameter_check_allowed_non_const.d \
./source/command/bcl_command_parameter_check_default.d \
./source/command/bcl_command_parameter_check_extension.d \
./source/command/bcl_command_parameter_check_extensions_file_existence.d \
./source/command/bcl_command_parameter_check_file_existence.d \
./source/command/bcl_command_parameter_check_file_in_search_path.d \
./source/command/bcl_command_parameter_check_interface.d \
./source/command/bcl_command_parameter_check_or.d \
./source/command/bcl_command_parameter_check_ranged.d \
./source/command/bcl_command_parameter_check_serializable.d 


# Each subdirectory must supply rules for building sources it contributes
source/command/%.o: ../source/command/%.cpp source/command/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



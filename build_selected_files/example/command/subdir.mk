################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/command/example_command.cpp \
../example/command/example_command_command.cpp \
../example/command/example_command_command_line_writer.cpp \
../example/command/example_command_flag_dynamic.cpp \
../example/command/example_command_flag_separator.cpp \
../example/command/example_command_flag_static.cpp \
../example/command/example_command_flag_static_and_dynamic.cpp \
../example/command/example_command_guesser.cpp \
../example/command/example_command_parameter.cpp \
../example/command/example_command_parameter_check_allowed.cpp \
../example/command/example_command_parameter_check_allowed_non_const.cpp \
../example/command/example_command_parameter_check_default.cpp \
../example/command/example_command_parameter_check_enumerate.cpp \
../example/command/example_command_parameter_check_extension.cpp \
../example/command/example_command_parameter_check_extensions_file_existence.cpp \
../example/command/example_command_parameter_check_file_existence.cpp \
../example/command/example_command_parameter_check_file_in_search_path.cpp \
../example/command/example_command_parameter_check_or.cpp \
../example/command/example_command_parameter_check_ranged.cpp \
../example/command/example_command_parameter_check_serializable.cpp 

OBJS += \
./example/command/example_command.o \
./example/command/example_command_command.o \
./example/command/example_command_command_line_writer.o \
./example/command/example_command_flag_dynamic.o \
./example/command/example_command_flag_separator.o \
./example/command/example_command_flag_static.o \
./example/command/example_command_flag_static_and_dynamic.o \
./example/command/example_command_guesser.o \
./example/command/example_command_parameter.o \
./example/command/example_command_parameter_check_allowed.o \
./example/command/example_command_parameter_check_allowed_non_const.o \
./example/command/example_command_parameter_check_default.o \
./example/command/example_command_parameter_check_enumerate.o \
./example/command/example_command_parameter_check_extension.o \
./example/command/example_command_parameter_check_extensions_file_existence.o \
./example/command/example_command_parameter_check_file_existence.o \
./example/command/example_command_parameter_check_file_in_search_path.o \
./example/command/example_command_parameter_check_or.o \
./example/command/example_command_parameter_check_ranged.o \
./example/command/example_command_parameter_check_serializable.o 

CPP_DEPS += \
./example/command/example_command.d \
./example/command/example_command_command.d \
./example/command/example_command_command_line_writer.d \
./example/command/example_command_flag_dynamic.d \
./example/command/example_command_flag_separator.d \
./example/command/example_command_flag_static.d \
./example/command/example_command_flag_static_and_dynamic.d \
./example/command/example_command_guesser.d \
./example/command/example_command_parameter.d \
./example/command/example_command_parameter_check_allowed.d \
./example/command/example_command_parameter_check_allowed_non_const.d \
./example/command/example_command_parameter_check_default.d \
./example/command/example_command_parameter_check_enumerate.d \
./example/command/example_command_parameter_check_extension.d \
./example/command/example_command_parameter_check_extensions_file_existence.d \
./example/command/example_command_parameter_check_file_existence.d \
./example/command/example_command_parameter_check_file_in_search_path.d \
./example/command/example_command_parameter_check_or.d \
./example/command/example_command_parameter_check_ranged.d \
./example/command/example_command_parameter_check_serializable.d 


# Each subdirectory must supply rules for building sources it contributes
example/command/%.o: ../example/command/%.cpp example/command/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



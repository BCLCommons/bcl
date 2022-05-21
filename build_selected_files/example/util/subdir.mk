################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/util/example_util.cpp \
../example/util/example_util_assert.cpp \
../example/util/example_util_call_stack.cpp \
../example/util/example_util_class_descriptor.cpp \
../example/util/example_util_color_gradient.cpp \
../example/util/example_util_colors.cpp \
../example/util/example_util_cpp_data_types.cpp \
../example/util/example_util_cpu_benchmark_whetstone.cpp \
../example/util/example_util_enum.cpp \
../example/util/example_util_format.cpp \
../example/util/example_util_implementation.cpp \
../example/util/example_util_logger_default.cpp \
../example/util/example_util_loggers.cpp \
../example/util/example_util_memory_usage.cpp \
../example/util/example_util_message.cpp \
../example/util/example_util_object_data_label.cpp \
../example/util/example_util_object_data_label_tokenizer.cpp \
../example/util/example_util_object_instances.cpp \
../example/util/example_util_object_interface.cpp \
../example/util/example_util_runtime_environment_default.cpp \
../example/util/example_util_runtime_environment_interface.cpp \
../example/util/example_util_runtime_environments.cpp \
../example/util/example_util_sh_ptr.cpp \
../example/util/example_util_sh_ptr_list.cpp \
../example/util/example_util_sh_ptr_vector.cpp \
../example/util/example_util_si_ptr.cpp \
../example/util/example_util_si_ptr_list.cpp \
../example/util/example_util_si_ptr_vector.cpp \
../example/util/example_util_stopwatch.cpp \
../example/util/example_util_string_functions.cpp \
../example/util/example_util_string_replacement.cpp \
../example/util/example_util_thunk_wrapper.cpp \
../example/util/example_util_time.cpp \
../example/util/example_util_undefined.cpp \
../example/util/example_util_wrapper.cpp \
../example/util/example_util_wrapper_base.cpp 

OBJS += \
./example/util/example_util.o \
./example/util/example_util_assert.o \
./example/util/example_util_call_stack.o \
./example/util/example_util_class_descriptor.o \
./example/util/example_util_color_gradient.o \
./example/util/example_util_colors.o \
./example/util/example_util_cpp_data_types.o \
./example/util/example_util_cpu_benchmark_whetstone.o \
./example/util/example_util_enum.o \
./example/util/example_util_format.o \
./example/util/example_util_implementation.o \
./example/util/example_util_logger_default.o \
./example/util/example_util_loggers.o \
./example/util/example_util_memory_usage.o \
./example/util/example_util_message.o \
./example/util/example_util_object_data_label.o \
./example/util/example_util_object_data_label_tokenizer.o \
./example/util/example_util_object_instances.o \
./example/util/example_util_object_interface.o \
./example/util/example_util_runtime_environment_default.o \
./example/util/example_util_runtime_environment_interface.o \
./example/util/example_util_runtime_environments.o \
./example/util/example_util_sh_ptr.o \
./example/util/example_util_sh_ptr_list.o \
./example/util/example_util_sh_ptr_vector.o \
./example/util/example_util_si_ptr.o \
./example/util/example_util_si_ptr_list.o \
./example/util/example_util_si_ptr_vector.o \
./example/util/example_util_stopwatch.o \
./example/util/example_util_string_functions.o \
./example/util/example_util_string_replacement.o \
./example/util/example_util_thunk_wrapper.o \
./example/util/example_util_time.o \
./example/util/example_util_undefined.o \
./example/util/example_util_wrapper.o \
./example/util/example_util_wrapper_base.o 

CPP_DEPS += \
./example/util/example_util.d \
./example/util/example_util_assert.d \
./example/util/example_util_call_stack.d \
./example/util/example_util_class_descriptor.d \
./example/util/example_util_color_gradient.d \
./example/util/example_util_colors.d \
./example/util/example_util_cpp_data_types.d \
./example/util/example_util_cpu_benchmark_whetstone.d \
./example/util/example_util_enum.d \
./example/util/example_util_format.d \
./example/util/example_util_implementation.d \
./example/util/example_util_logger_default.d \
./example/util/example_util_loggers.d \
./example/util/example_util_memory_usage.d \
./example/util/example_util_message.d \
./example/util/example_util_object_data_label.d \
./example/util/example_util_object_data_label_tokenizer.d \
./example/util/example_util_object_instances.d \
./example/util/example_util_object_interface.d \
./example/util/example_util_runtime_environment_default.d \
./example/util/example_util_runtime_environment_interface.d \
./example/util/example_util_runtime_environments.d \
./example/util/example_util_sh_ptr.d \
./example/util/example_util_sh_ptr_list.d \
./example/util/example_util_sh_ptr_vector.d \
./example/util/example_util_si_ptr.d \
./example/util/example_util_si_ptr_list.d \
./example/util/example_util_si_ptr_vector.d \
./example/util/example_util_stopwatch.d \
./example/util/example_util_string_functions.d \
./example/util/example_util_string_replacement.d \
./example/util/example_util_thunk_wrapper.d \
./example/util/example_util_time.d \
./example/util/example_util_undefined.d \
./example/util/example_util_wrapper.d \
./example/util/example_util_wrapper_base.d 


# Each subdirectory must supply rules for building sources it contributes
example/util/%.o: ../example/util/%.cpp example/util/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



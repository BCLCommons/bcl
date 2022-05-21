################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/sched/example_sched.cpp \
../example/sched/example_sched_function_job_with_data.cpp \
../example/sched/example_sched_jobs_with_data.cpp \
../example/sched/example_sched_mutex.cpp \
../example/sched/example_sched_serial_scheduler.cpp \
../example/sched/example_sched_sum_function.cpp 

OBJS += \
./example/sched/example_sched.o \
./example/sched/example_sched_function_job_with_data.o \
./example/sched/example_sched_jobs_with_data.o \
./example/sched/example_sched_mutex.o \
./example/sched/example_sched_serial_scheduler.o \
./example/sched/example_sched_sum_function.o 

CPP_DEPS += \
./example/sched/example_sched.d \
./example/sched/example_sched_function_job_with_data.d \
./example/sched/example_sched_jobs_with_data.d \
./example/sched/example_sched_mutex.d \
./example/sched/example_sched_serial_scheduler.d \
./example/sched/example_sched_sum_function.d 


# Each subdirectory must supply rules for building sources it contributes
example/sched/%.o: ../example/sched/%.cpp example/sched/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



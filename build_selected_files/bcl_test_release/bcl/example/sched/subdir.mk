################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/sched/example_sched.cpp \
../bcl_test_release/bcl/example/sched/example_sched_function_job_with_data.cpp \
../bcl_test_release/bcl/example/sched/example_sched_jobs_with_data.cpp \
../bcl_test_release/bcl/example/sched/example_sched_mutex.cpp \
../bcl_test_release/bcl/example/sched/example_sched_serial_scheduler.cpp \
../bcl_test_release/bcl/example/sched/example_sched_sum_function.cpp 

OBJS += \
./bcl_test_release/bcl/example/sched/example_sched.o \
./bcl_test_release/bcl/example/sched/example_sched_function_job_with_data.o \
./bcl_test_release/bcl/example/sched/example_sched_jobs_with_data.o \
./bcl_test_release/bcl/example/sched/example_sched_mutex.o \
./bcl_test_release/bcl/example/sched/example_sched_serial_scheduler.o \
./bcl_test_release/bcl/example/sched/example_sched_sum_function.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/sched/example_sched.d \
./bcl_test_release/bcl/example/sched/example_sched_function_job_with_data.d \
./bcl_test_release/bcl/example/sched/example_sched_jobs_with_data.d \
./bcl_test_release/bcl/example/sched/example_sched_mutex.d \
./bcl_test_release/bcl/example/sched/example_sched_serial_scheduler.d \
./bcl_test_release/bcl/example/sched/example_sched_sum_function.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/sched/%.o: ../bcl_test_release/bcl/example/sched/%.cpp bcl_test_release/bcl/example/sched/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



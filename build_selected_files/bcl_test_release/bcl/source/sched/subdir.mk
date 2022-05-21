################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/sched/bcl_sched.cpp \
../bcl_test_release/bcl/source/sched/bcl_sched_mutex.cpp \
../bcl_test_release/bcl/source/sched/bcl_sched_scheduler_interface.cpp \
../bcl_test_release/bcl/source/sched/bcl_sched_schedulers.cpp 

OBJS += \
./bcl_test_release/bcl/source/sched/bcl_sched.o \
./bcl_test_release/bcl/source/sched/bcl_sched_mutex.o \
./bcl_test_release/bcl/source/sched/bcl_sched_scheduler_interface.o \
./bcl_test_release/bcl/source/sched/bcl_sched_schedulers.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/sched/bcl_sched.d \
./bcl_test_release/bcl/source/sched/bcl_sched_mutex.d \
./bcl_test_release/bcl/source/sched/bcl_sched_scheduler_interface.d \
./bcl_test_release/bcl/source/sched/bcl_sched_schedulers.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/sched/%.o: ../bcl_test_release/bcl/source/sched/%.cpp bcl_test_release/bcl/source/sched/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



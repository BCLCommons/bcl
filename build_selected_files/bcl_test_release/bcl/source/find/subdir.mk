################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/find/bcl_find.cpp \
../bcl_test_release/bcl/source/find/bcl_find_pick_body_extent.cpp \
../bcl_test_release/bcl/source/find/bcl_find_pick_body_random.cpp \
../bcl_test_release/bcl/source/find/bcl_find_template_instantiations.cpp 

OBJS += \
./bcl_test_release/bcl/source/find/bcl_find.o \
./bcl_test_release/bcl/source/find/bcl_find_pick_body_extent.o \
./bcl_test_release/bcl/source/find/bcl_find_pick_body_random.o \
./bcl_test_release/bcl/source/find/bcl_find_template_instantiations.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/find/bcl_find.d \
./bcl_test_release/bcl/source/find/bcl_find_pick_body_extent.d \
./bcl_test_release/bcl/source/find/bcl_find_pick_body_random.d \
./bcl_test_release/bcl/source/find/bcl_find_template_instantiations.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/find/%.o: ../bcl_test_release/bcl/source/find/%.cpp bcl_test_release/bcl/source/find/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



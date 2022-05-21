################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/apps/internal/biology/example_app_analyze_loops.cpp \
../bcl_test_release/bcl/example/apps/internal/biology/example_app_loop_template.cpp \
../bcl_test_release/bcl/example/apps/internal/biology/example_app_scatterplot_from_tables.cpp 

OBJS += \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_analyze_loops.o \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_loop_template.o \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_scatterplot_from_tables.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_analyze_loops.d \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_loop_template.d \
./bcl_test_release/bcl/example/apps/internal/biology/example_app_scatterplot_from_tables.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/apps/internal/biology/%.o: ../bcl_test_release/bcl/example/apps/internal/biology/%.cpp bcl_test_release/bcl/example/apps/internal/biology/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



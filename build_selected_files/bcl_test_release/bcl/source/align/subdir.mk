################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/align/bcl_align.cpp \
../bcl_test_release/bcl/source/align/bcl_align_handler_classes.cpp \
../bcl_test_release/bcl/source/align/bcl_align_multiple_aligner_classes.cpp \
../bcl_test_release/bcl/source/align/bcl_align_pairwise_aligner_classes.cpp 

OBJS += \
./bcl_test_release/bcl/source/align/bcl_align.o \
./bcl_test_release/bcl/source/align/bcl_align_handler_classes.o \
./bcl_test_release/bcl/source/align/bcl_align_multiple_aligner_classes.o \
./bcl_test_release/bcl/source/align/bcl_align_pairwise_aligner_classes.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/align/bcl_align.d \
./bcl_test_release/bcl/source/align/bcl_align_handler_classes.d \
./bcl_test_release/bcl/source/align/bcl_align_multiple_aligner_classes.d \
./bcl_test_release/bcl/source/align/bcl_align_pairwise_aligner_classes.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/align/%.o: ../bcl_test_release/bcl/source/align/%.cpp bcl_test_release/bcl/source/align/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



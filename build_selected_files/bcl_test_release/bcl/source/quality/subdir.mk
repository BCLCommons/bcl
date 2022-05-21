################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/quality/bcl_quality.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_average.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_const_measure.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_dme.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_dmf.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_gdt.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_lcs.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_maxsub.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_measures.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_rmsd.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_rmsd_preprocessor.cpp \
../bcl_test_release/bcl/source/quality/bcl_quality_superimpose_measures.cpp 

OBJS += \
./bcl_test_release/bcl/source/quality/bcl_quality.o \
./bcl_test_release/bcl/source/quality/bcl_quality_average.o \
./bcl_test_release/bcl/source/quality/bcl_quality_const_measure.o \
./bcl_test_release/bcl/source/quality/bcl_quality_dme.o \
./bcl_test_release/bcl/source/quality/bcl_quality_dmf.o \
./bcl_test_release/bcl/source/quality/bcl_quality_gdt.o \
./bcl_test_release/bcl/source/quality/bcl_quality_lcs.o \
./bcl_test_release/bcl/source/quality/bcl_quality_maxsub.o \
./bcl_test_release/bcl/source/quality/bcl_quality_measures.o \
./bcl_test_release/bcl/source/quality/bcl_quality_rmsd.o \
./bcl_test_release/bcl/source/quality/bcl_quality_rmsd_preprocessor.o \
./bcl_test_release/bcl/source/quality/bcl_quality_superimpose_measures.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/quality/bcl_quality.d \
./bcl_test_release/bcl/source/quality/bcl_quality_average.d \
./bcl_test_release/bcl/source/quality/bcl_quality_const_measure.d \
./bcl_test_release/bcl/source/quality/bcl_quality_dme.d \
./bcl_test_release/bcl/source/quality/bcl_quality_dmf.d \
./bcl_test_release/bcl/source/quality/bcl_quality_gdt.d \
./bcl_test_release/bcl/source/quality/bcl_quality_lcs.d \
./bcl_test_release/bcl/source/quality/bcl_quality_maxsub.d \
./bcl_test_release/bcl/source/quality/bcl_quality_measures.d \
./bcl_test_release/bcl/source/quality/bcl_quality_rmsd.d \
./bcl_test_release/bcl/source/quality/bcl_quality_rmsd_preprocessor.d \
./bcl_test_release/bcl/source/quality/bcl_quality_superimpose_measures.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/quality/%.o: ../bcl_test_release/bcl/source/quality/%.cpp bcl_test_release/bcl/source/quality/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



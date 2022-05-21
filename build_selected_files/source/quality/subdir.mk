################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/quality/bcl_quality.cpp \
../source/quality/bcl_quality_average.cpp \
../source/quality/bcl_quality_const_measure.cpp \
../source/quality/bcl_quality_dme.cpp \
../source/quality/bcl_quality_dmf.cpp \
../source/quality/bcl_quality_gdt.cpp \
../source/quality/bcl_quality_lcs.cpp \
../source/quality/bcl_quality_maxsub.cpp \
../source/quality/bcl_quality_measures.cpp \
../source/quality/bcl_quality_rmsd.cpp \
../source/quality/bcl_quality_rmsd_preprocessor.cpp \
../source/quality/bcl_quality_superimpose_measures.cpp 

OBJS += \
./source/quality/bcl_quality.o \
./source/quality/bcl_quality_average.o \
./source/quality/bcl_quality_const_measure.o \
./source/quality/bcl_quality_dme.o \
./source/quality/bcl_quality_dmf.o \
./source/quality/bcl_quality_gdt.o \
./source/quality/bcl_quality_lcs.o \
./source/quality/bcl_quality_maxsub.o \
./source/quality/bcl_quality_measures.o \
./source/quality/bcl_quality_rmsd.o \
./source/quality/bcl_quality_rmsd_preprocessor.o \
./source/quality/bcl_quality_superimpose_measures.o 

CPP_DEPS += \
./source/quality/bcl_quality.d \
./source/quality/bcl_quality_average.d \
./source/quality/bcl_quality_const_measure.d \
./source/quality/bcl_quality_dme.d \
./source/quality/bcl_quality_dmf.d \
./source/quality/bcl_quality_gdt.d \
./source/quality/bcl_quality_lcs.d \
./source/quality/bcl_quality_maxsub.d \
./source/quality/bcl_quality_measures.d \
./source/quality/bcl_quality_rmsd.d \
./source/quality/bcl_quality_rmsd_preprocessor.d \
./source/quality/bcl_quality_superimpose_measures.d 


# Each subdirectory must supply rules for building sources it contributes
source/quality/%.o: ../source/quality/%.cpp source/quality/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/quality/example_quality.cpp \
../example/quality/example_quality_const_measure.cpp \
../example/quality/example_quality_dme.cpp \
../example/quality/example_quality_dmf.cpp \
../example/quality/example_quality_gdt.cpp \
../example/quality/example_quality_lcs.cpp \
../example/quality/example_quality_maxsub.cpp \
../example/quality/example_quality_measures.cpp \
../example/quality/example_quality_rmsd.cpp \
../example/quality/example_quality_superimpose_measures.cpp 

OBJS += \
./example/quality/example_quality.o \
./example/quality/example_quality_const_measure.o \
./example/quality/example_quality_dme.o \
./example/quality/example_quality_dmf.o \
./example/quality/example_quality_gdt.o \
./example/quality/example_quality_lcs.o \
./example/quality/example_quality_maxsub.o \
./example/quality/example_quality_measures.o \
./example/quality/example_quality_rmsd.o \
./example/quality/example_quality_superimpose_measures.o 

CPP_DEPS += \
./example/quality/example_quality.d \
./example/quality/example_quality_const_measure.d \
./example/quality/example_quality_dme.d \
./example/quality/example_quality_dmf.d \
./example/quality/example_quality_gdt.d \
./example/quality/example_quality_lcs.d \
./example/quality/example_quality_maxsub.d \
./example/quality/example_quality_measures.d \
./example/quality/example_quality_rmsd.d \
./example/quality/example_quality_superimpose_measures.d 


# Each subdirectory must supply rules for building sources it contributes
example/quality/%.o: ../example/quality/%.cpp example/quality/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



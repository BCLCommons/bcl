################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/sspred/example_sspred.cpp \
../example/sspred/example_sspred_b2tmpred.cpp \
../example/sspred/example_sspred_boctopus.cpp \
../example/sspred/example_sspred_conpred.cpp \
../example/sspred/example_sspred_dssp.cpp \
../example/sspred/example_sspred_dssp_stride.cpp \
../example/sspred/example_sspred_jufo.cpp \
../example/sspred/example_sspred_jufo9d.cpp \
../example/sspred/example_sspred_kaksi.cpp \
../example/sspred/example_sspred_masp.cpp \
../example/sspred/example_sspred_method_handler.cpp \
../example/sspred/example_sspred_methods.cpp \
../example/sspred/example_sspred_octopus.cpp \
../example/sspred/example_sspred_palsse.cpp \
../example/sspred/example_sspred_partifold.cpp \
../example/sspred/example_sspred_pdb.cpp \
../example/sspred/example_sspred_profphd.cpp \
../example/sspred/example_sspred_proftmb.cpp \
../example/sspred/example_sspred_psipred.cpp \
../example/sspred/example_sspred_sam.cpp \
../example/sspred/example_sspred_sse_factory_highest.cpp \
../example/sspred/example_sspred_sse_factory_threshold.cpp \
../example/sspred/example_sspred_stride.cpp \
../example/sspred/example_sspred_talos.cpp \
../example/sspred/example_sspred_tmbetanet.cpp \
../example/sspred/example_sspred_tmhmm.cpp \
../example/sspred/example_sspred_tmmod.cpp 

OBJS += \
./example/sspred/example_sspred.o \
./example/sspred/example_sspred_b2tmpred.o \
./example/sspred/example_sspred_boctopus.o \
./example/sspred/example_sspred_conpred.o \
./example/sspred/example_sspred_dssp.o \
./example/sspred/example_sspred_dssp_stride.o \
./example/sspred/example_sspred_jufo.o \
./example/sspred/example_sspred_jufo9d.o \
./example/sspred/example_sspred_kaksi.o \
./example/sspred/example_sspred_masp.o \
./example/sspred/example_sspred_method_handler.o \
./example/sspred/example_sspred_methods.o \
./example/sspred/example_sspred_octopus.o \
./example/sspred/example_sspred_palsse.o \
./example/sspred/example_sspred_partifold.o \
./example/sspred/example_sspred_pdb.o \
./example/sspred/example_sspred_profphd.o \
./example/sspred/example_sspred_proftmb.o \
./example/sspred/example_sspred_psipred.o \
./example/sspred/example_sspred_sam.o \
./example/sspred/example_sspred_sse_factory_highest.o \
./example/sspred/example_sspred_sse_factory_threshold.o \
./example/sspred/example_sspred_stride.o \
./example/sspred/example_sspred_talos.o \
./example/sspred/example_sspred_tmbetanet.o \
./example/sspred/example_sspred_tmhmm.o \
./example/sspred/example_sspred_tmmod.o 

CPP_DEPS += \
./example/sspred/example_sspred.d \
./example/sspred/example_sspred_b2tmpred.d \
./example/sspred/example_sspred_boctopus.d \
./example/sspred/example_sspred_conpred.d \
./example/sspred/example_sspred_dssp.d \
./example/sspred/example_sspred_dssp_stride.d \
./example/sspred/example_sspred_jufo.d \
./example/sspred/example_sspred_jufo9d.d \
./example/sspred/example_sspred_kaksi.d \
./example/sspred/example_sspred_masp.d \
./example/sspred/example_sspred_method_handler.d \
./example/sspred/example_sspred_methods.d \
./example/sspred/example_sspred_octopus.d \
./example/sspred/example_sspred_palsse.d \
./example/sspred/example_sspred_partifold.d \
./example/sspred/example_sspred_pdb.d \
./example/sspred/example_sspred_profphd.d \
./example/sspred/example_sspred_proftmb.d \
./example/sspred/example_sspred_psipred.d \
./example/sspred/example_sspred_sam.d \
./example/sspred/example_sspred_sse_factory_highest.d \
./example/sspred/example_sspred_sse_factory_threshold.d \
./example/sspred/example_sspred_stride.d \
./example/sspred/example_sspred_talos.d \
./example/sspred/example_sspred_tmbetanet.d \
./example/sspred/example_sspred_tmhmm.d \
./example/sspred/example_sspred_tmmod.d 


# Each subdirectory must supply rules for building sources it contributes
example/sspred/%.o: ../example/sspred/%.cpp example/sspred/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



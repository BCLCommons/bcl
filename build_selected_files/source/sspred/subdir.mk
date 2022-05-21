################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/sspred/bcl_sspred.cpp \
../source/sspred/bcl_sspred_b2tmpred.cpp \
../source/sspred/bcl_sspred_boctopus.cpp \
../source/sspred/bcl_sspred_ci_phi_psi.cpp \
../source/sspred/bcl_sspred_conpred.cpp \
../source/sspred/bcl_sspred_dssp.cpp \
../source/sspred/bcl_sspred_dssp_stride.cpp \
../source/sspred/bcl_sspred_jufo.cpp \
../source/sspred/bcl_sspred_jufo9d.cpp \
../source/sspred/bcl_sspred_jufo_ann.cpp \
../source/sspred/bcl_sspred_kaksi.cpp \
../source/sspred/bcl_sspred_mahssmi.cpp \
../source/sspred/bcl_sspred_masp.cpp \
../source/sspred/bcl_sspred_method_handler.cpp \
../source/sspred/bcl_sspred_method_interface.cpp \
../source/sspred/bcl_sspred_methods.cpp \
../source/sspred/bcl_sspred_octopus.cpp \
../source/sspred/bcl_sspred_palsse.cpp \
../source/sspred/bcl_sspred_partifold.cpp \
../source/sspred/bcl_sspred_pdb.cpp \
../source/sspred/bcl_sspred_profphd.cpp \
../source/sspred/bcl_sspred_proftmb.cpp \
../source/sspred/bcl_sspred_psipred.cpp \
../source/sspred/bcl_sspred_sam.cpp \
../source/sspred/bcl_sspred_sse_factory_highest.cpp \
../source/sspred/bcl_sspred_sse_factory_threshold.cpp \
../source/sspred/bcl_sspred_stride.cpp \
../source/sspred/bcl_sspred_talos.cpp \
../source/sspred/bcl_sspred_tmbetanet.cpp \
../source/sspred/bcl_sspred_tmhmm.cpp \
../source/sspred/bcl_sspred_tmmod.cpp 

OBJS += \
./source/sspred/bcl_sspred.o \
./source/sspred/bcl_sspred_b2tmpred.o \
./source/sspred/bcl_sspred_boctopus.o \
./source/sspred/bcl_sspred_ci_phi_psi.o \
./source/sspred/bcl_sspred_conpred.o \
./source/sspred/bcl_sspred_dssp.o \
./source/sspred/bcl_sspred_dssp_stride.o \
./source/sspred/bcl_sspred_jufo.o \
./source/sspred/bcl_sspred_jufo9d.o \
./source/sspred/bcl_sspred_jufo_ann.o \
./source/sspred/bcl_sspred_kaksi.o \
./source/sspred/bcl_sspred_mahssmi.o \
./source/sspred/bcl_sspred_masp.o \
./source/sspred/bcl_sspred_method_handler.o \
./source/sspred/bcl_sspred_method_interface.o \
./source/sspred/bcl_sspred_methods.o \
./source/sspred/bcl_sspred_octopus.o \
./source/sspred/bcl_sspred_palsse.o \
./source/sspred/bcl_sspred_partifold.o \
./source/sspred/bcl_sspred_pdb.o \
./source/sspred/bcl_sspred_profphd.o \
./source/sspred/bcl_sspred_proftmb.o \
./source/sspred/bcl_sspred_psipred.o \
./source/sspred/bcl_sspred_sam.o \
./source/sspred/bcl_sspred_sse_factory_highest.o \
./source/sspred/bcl_sspred_sse_factory_threshold.o \
./source/sspred/bcl_sspred_stride.o \
./source/sspred/bcl_sspred_talos.o \
./source/sspred/bcl_sspred_tmbetanet.o \
./source/sspred/bcl_sspred_tmhmm.o \
./source/sspred/bcl_sspred_tmmod.o 

CPP_DEPS += \
./source/sspred/bcl_sspred.d \
./source/sspred/bcl_sspred_b2tmpred.d \
./source/sspred/bcl_sspred_boctopus.d \
./source/sspred/bcl_sspred_ci_phi_psi.d \
./source/sspred/bcl_sspred_conpred.d \
./source/sspred/bcl_sspred_dssp.d \
./source/sspred/bcl_sspred_dssp_stride.d \
./source/sspred/bcl_sspred_jufo.d \
./source/sspred/bcl_sspred_jufo9d.d \
./source/sspred/bcl_sspred_jufo_ann.d \
./source/sspred/bcl_sspred_kaksi.d \
./source/sspred/bcl_sspred_mahssmi.d \
./source/sspred/bcl_sspred_masp.d \
./source/sspred/bcl_sspred_method_handler.d \
./source/sspred/bcl_sspred_method_interface.d \
./source/sspred/bcl_sspred_methods.d \
./source/sspred/bcl_sspred_octopus.d \
./source/sspred/bcl_sspred_palsse.d \
./source/sspred/bcl_sspred_partifold.d \
./source/sspred/bcl_sspred_pdb.d \
./source/sspred/bcl_sspred_profphd.d \
./source/sspred/bcl_sspred_proftmb.d \
./source/sspred/bcl_sspred_psipred.d \
./source/sspred/bcl_sspred_sam.d \
./source/sspred/bcl_sspred_sse_factory_highest.d \
./source/sspred/bcl_sspred_sse_factory_threshold.d \
./source/sspred/bcl_sspred_stride.d \
./source/sspred/bcl_sspred_talos.d \
./source/sspred/bcl_sspred_tmbetanet.d \
./source/sspred/bcl_sspred_tmhmm.d \
./source/sspred/bcl_sspred_tmmod.d 


# Each subdirectory must supply rules for building sources it contributes
source/sspred/%.o: ../source/sspred/%.cpp source/sspred/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/release/bcl_app_alignment.cpp \
../apps/release/bcl_app_cluster.cpp \
../apps/release/bcl_app_contact_prediction.cpp \
../apps/release/bcl_app_create_sse_pool.cpp \
../apps/release/bcl_app_fit_in_density.cpp \
../apps/release/bcl_app_fit_in_density_minimize.cpp \
../apps/release/bcl_app_fold.cpp \
../apps/release/bcl_app_fusion_protein.cpp \
../apps/release/bcl_app_jufo.cpp \
../apps/release/bcl_app_optimize.cpp \
../apps/release/bcl_app_pdb_compare.cpp \
../apps/release/bcl_app_pdb_convert.cpp \
../apps/release/bcl_app_protein_docking.cpp \
../apps/release/bcl_app_score_protein.cpp 

OBJS += \
./apps/release/bcl_app_alignment.o \
./apps/release/bcl_app_cluster.o \
./apps/release/bcl_app_contact_prediction.o \
./apps/release/bcl_app_create_sse_pool.o \
./apps/release/bcl_app_fit_in_density.o \
./apps/release/bcl_app_fit_in_density_minimize.o \
./apps/release/bcl_app_fold.o \
./apps/release/bcl_app_fusion_protein.o \
./apps/release/bcl_app_jufo.o \
./apps/release/bcl_app_optimize.o \
./apps/release/bcl_app_pdb_compare.o \
./apps/release/bcl_app_pdb_convert.o \
./apps/release/bcl_app_protein_docking.o \
./apps/release/bcl_app_score_protein.o 

CPP_DEPS += \
./apps/release/bcl_app_alignment.d \
./apps/release/bcl_app_cluster.d \
./apps/release/bcl_app_contact_prediction.d \
./apps/release/bcl_app_create_sse_pool.d \
./apps/release/bcl_app_fit_in_density.d \
./apps/release/bcl_app_fit_in_density_minimize.d \
./apps/release/bcl_app_fold.d \
./apps/release/bcl_app_fusion_protein.d \
./apps/release/bcl_app_jufo.d \
./apps/release/bcl_app_optimize.d \
./apps/release/bcl_app_pdb_compare.d \
./apps/release/bcl_app_pdb_convert.d \
./apps/release/bcl_app_protein_docking.d \
./apps/release/bcl_app_score_protein.d 


# Each subdirectory must supply rules for building sources it contributes
apps/release/%.o: ../apps/release/%.cpp apps/release/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



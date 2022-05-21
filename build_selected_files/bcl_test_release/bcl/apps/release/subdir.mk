################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/apps/release/bcl_app_alignment.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_cluster.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_contact_prediction.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_create_sse_pool.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_fit_in_density.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_fit_in_density_minimize.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_fold.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_fusion_protein.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_jufo.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_optimize.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_pdb_compare.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_pdb_convert.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_protein_docking.cpp \
../bcl_test_release/bcl/apps/release/bcl_app_score_protein.cpp 

OBJS += \
./bcl_test_release/bcl/apps/release/bcl_app_alignment.o \
./bcl_test_release/bcl/apps/release/bcl_app_cluster.o \
./bcl_test_release/bcl/apps/release/bcl_app_contact_prediction.o \
./bcl_test_release/bcl/apps/release/bcl_app_create_sse_pool.o \
./bcl_test_release/bcl/apps/release/bcl_app_fit_in_density.o \
./bcl_test_release/bcl/apps/release/bcl_app_fit_in_density_minimize.o \
./bcl_test_release/bcl/apps/release/bcl_app_fold.o \
./bcl_test_release/bcl/apps/release/bcl_app_fusion_protein.o \
./bcl_test_release/bcl/apps/release/bcl_app_jufo.o \
./bcl_test_release/bcl/apps/release/bcl_app_optimize.o \
./bcl_test_release/bcl/apps/release/bcl_app_pdb_compare.o \
./bcl_test_release/bcl/apps/release/bcl_app_pdb_convert.o \
./bcl_test_release/bcl/apps/release/bcl_app_protein_docking.o \
./bcl_test_release/bcl/apps/release/bcl_app_score_protein.o 

CPP_DEPS += \
./bcl_test_release/bcl/apps/release/bcl_app_alignment.d \
./bcl_test_release/bcl/apps/release/bcl_app_cluster.d \
./bcl_test_release/bcl/apps/release/bcl_app_contact_prediction.d \
./bcl_test_release/bcl/apps/release/bcl_app_create_sse_pool.d \
./bcl_test_release/bcl/apps/release/bcl_app_fit_in_density.d \
./bcl_test_release/bcl/apps/release/bcl_app_fit_in_density_minimize.d \
./bcl_test_release/bcl/apps/release/bcl_app_fold.d \
./bcl_test_release/bcl/apps/release/bcl_app_fusion_protein.d \
./bcl_test_release/bcl/apps/release/bcl_app_jufo.d \
./bcl_test_release/bcl/apps/release/bcl_app_optimize.d \
./bcl_test_release/bcl/apps/release/bcl_app_pdb_compare.d \
./bcl_test_release/bcl/apps/release/bcl_app_pdb_convert.d \
./bcl_test_release/bcl/apps/release/bcl_app_protein_docking.d \
./bcl_test_release/bcl/apps/release/bcl_app_score_protein.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/apps/release/%.o: ../bcl_test_release/bcl/apps/release/%.cpp bcl_test_release/bcl/apps/release/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



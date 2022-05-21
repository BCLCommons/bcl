################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/apps/release/example_app_alignment.cpp \
../example/apps/release/example_app_cluster.cpp \
../example/apps/release/example_app_contact_prediction.cpp \
../example/apps/release/example_app_create_sse_pool.cpp \
../example/apps/release/example_app_fit_in_density.cpp \
../example/apps/release/example_app_fit_in_density_minimize.cpp \
../example/apps/release/example_app_fold.cpp \
../example/apps/release/example_app_fusion_protein.cpp \
../example/apps/release/example_app_jufo.cpp \
../example/apps/release/example_app_optimize.cpp \
../example/apps/release/example_app_pdb_compare.cpp \
../example/apps/release/example_app_pdb_convert.cpp \
../example/apps/release/example_app_score_protein.cpp 

OBJS += \
./example/apps/release/example_app_alignment.o \
./example/apps/release/example_app_cluster.o \
./example/apps/release/example_app_contact_prediction.o \
./example/apps/release/example_app_create_sse_pool.o \
./example/apps/release/example_app_fit_in_density.o \
./example/apps/release/example_app_fit_in_density_minimize.o \
./example/apps/release/example_app_fold.o \
./example/apps/release/example_app_fusion_protein.o \
./example/apps/release/example_app_jufo.o \
./example/apps/release/example_app_optimize.o \
./example/apps/release/example_app_pdb_compare.o \
./example/apps/release/example_app_pdb_convert.o \
./example/apps/release/example_app_score_protein.o 

CPP_DEPS += \
./example/apps/release/example_app_alignment.d \
./example/apps/release/example_app_cluster.d \
./example/apps/release/example_app_contact_prediction.d \
./example/apps/release/example_app_create_sse_pool.d \
./example/apps/release/example_app_fit_in_density.d \
./example/apps/release/example_app_fit_in_density_minimize.d \
./example/apps/release/example_app_fold.d \
./example/apps/release/example_app_fusion_protein.d \
./example/apps/release/example_app_jufo.d \
./example/apps/release/example_app_optimize.d \
./example/apps/release/example_app_pdb_compare.d \
./example/apps/release/example_app_pdb_convert.d \
./example/apps/release/example_app_score_protein.d 


# Each subdirectory must supply rules for building sources it contributes
example/apps/release/%.o: ../example/apps/release/%.cpp example/apps/release/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



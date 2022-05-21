################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/scorestat/example_scorestat_loop_closure.cpp \
../example/scorestat/example_scorestat_loop_distance.cpp \
../example/scorestat/example_scorestat_neighbor_count.cpp \
../example/scorestat/example_scorestat_neighbor_vector.cpp \
../example/scorestat/example_scorestat_ols.cpp \
../example/scorestat/example_scorestat_phipsi.cpp \
../example/scorestat/example_scorestat_protein_model_packing.cpp \
../example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.cpp \
../example/scorestat/example_scorestat_radius_of_gyration.cpp 

OBJS += \
./example/scorestat/example_scorestat_loop_closure.o \
./example/scorestat/example_scorestat_loop_distance.o \
./example/scorestat/example_scorestat_neighbor_count.o \
./example/scorestat/example_scorestat_neighbor_vector.o \
./example/scorestat/example_scorestat_ols.o \
./example/scorestat/example_scorestat_phipsi.o \
./example/scorestat/example_scorestat_protein_model_packing.o \
./example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.o \
./example/scorestat/example_scorestat_radius_of_gyration.o 

CPP_DEPS += \
./example/scorestat/example_scorestat_loop_closure.d \
./example/scorestat/example_scorestat_loop_distance.d \
./example/scorestat/example_scorestat_neighbor_count.d \
./example/scorestat/example_scorestat_neighbor_vector.d \
./example/scorestat/example_scorestat_ols.d \
./example/scorestat/example_scorestat_phipsi.d \
./example/scorestat/example_scorestat_protein_model_packing.d \
./example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.d \
./example/scorestat/example_scorestat_radius_of_gyration.d 


# Each subdirectory must supply rules for building sources it contributes
example/scorestat/%.o: ../example/scorestat/%.cpp example/scorestat/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



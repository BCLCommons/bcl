################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/scorestat/example_scorestat_loop_closure.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_loop_distance.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_count.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_vector.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_ols.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_phipsi.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_packing.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.cpp \
../bcl_test_release/bcl/example/scorestat/example_scorestat_radius_of_gyration.cpp 

OBJS += \
./bcl_test_release/bcl/example/scorestat/example_scorestat_loop_closure.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_loop_distance.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_count.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_vector.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_ols.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_phipsi.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_packing.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.o \
./bcl_test_release/bcl/example/scorestat/example_scorestat_radius_of_gyration.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/scorestat/example_scorestat_loop_closure.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_loop_distance.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_count.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_neighbor_vector.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_ols.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_phipsi.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_packing.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_protein_model_sse_triplet_chirality.d \
./bcl_test_release/bcl/example/scorestat/example_scorestat_radius_of_gyration.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/scorestat/%.o: ../bcl_test_release/bcl/example/scorestat/%.cpp bcl_test_release/bcl/example/scorestat/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



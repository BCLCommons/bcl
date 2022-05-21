################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/source/cluster/bcl_cluster.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_input_classes.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_input_features.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_input_pairwise_list.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_input_table.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_linkage_classes.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_node_description_from_file.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_output_classes.cpp \
../bcl_test_release/bcl/source/cluster/bcl_cluster_output_pymol_label_string.cpp 

OBJS += \
./bcl_test_release/bcl/source/cluster/bcl_cluster.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_classes.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_features.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_pairwise_list.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_table.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_linkage_classes.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_node_description_from_file.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_output_classes.o \
./bcl_test_release/bcl/source/cluster/bcl_cluster_output_pymol_label_string.o 

CPP_DEPS += \
./bcl_test_release/bcl/source/cluster/bcl_cluster.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_classes.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_features.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_pairwise_list.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_input_table.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_linkage_classes.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_node_description_from_file.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_output_classes.d \
./bcl_test_release/bcl/source/cluster/bcl_cluster_output_pymol_label_string.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/source/cluster/%.o: ../bcl_test_release/bcl/source/cluster/%.cpp bcl_test_release/bcl/source/cluster/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/cluster/bcl_cluster.cpp \
../source/cluster/bcl_cluster_input_classes.cpp \
../source/cluster/bcl_cluster_input_features.cpp \
../source/cluster/bcl_cluster_input_pairwise_list.cpp \
../source/cluster/bcl_cluster_input_table.cpp \
../source/cluster/bcl_cluster_linkage_classes.cpp \
../source/cluster/bcl_cluster_node_description_from_file.cpp \
../source/cluster/bcl_cluster_output_classes.cpp \
../source/cluster/bcl_cluster_output_pymol_label_string.cpp 

OBJS += \
./source/cluster/bcl_cluster.o \
./source/cluster/bcl_cluster_input_classes.o \
./source/cluster/bcl_cluster_input_features.o \
./source/cluster/bcl_cluster_input_pairwise_list.o \
./source/cluster/bcl_cluster_input_table.o \
./source/cluster/bcl_cluster_linkage_classes.o \
./source/cluster/bcl_cluster_node_description_from_file.o \
./source/cluster/bcl_cluster_output_classes.o \
./source/cluster/bcl_cluster_output_pymol_label_string.o 

CPP_DEPS += \
./source/cluster/bcl_cluster.d \
./source/cluster/bcl_cluster_input_classes.d \
./source/cluster/bcl_cluster_input_features.d \
./source/cluster/bcl_cluster_input_pairwise_list.d \
./source/cluster/bcl_cluster_input_table.d \
./source/cluster/bcl_cluster_linkage_classes.d \
./source/cluster/bcl_cluster_node_description_from_file.d \
./source/cluster/bcl_cluster_output_classes.d \
./source/cluster/bcl_cluster_output_pymol_label_string.d 


# Each subdirectory must supply rules for building sources it contributes
source/cluster/%.o: ../source/cluster/%.cpp source/cluster/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



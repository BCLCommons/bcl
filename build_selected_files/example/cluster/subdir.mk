################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/cluster/example_cluster.cpp \
../example/cluster/example_cluster_dendrogram.cpp \
../example/cluster/example_cluster_distances_euclidean.cpp \
../example/cluster/example_cluster_distances_stored.cpp \
../example/cluster/example_cluster_input_classes.cpp \
../example/cluster/example_cluster_input_pairwise_list.cpp \
../example/cluster/example_cluster_input_table.cpp \
../example/cluster/example_cluster_linkage_average.cpp \
../example/cluster/example_cluster_linkage_classes.cpp \
../example/cluster/example_cluster_linkage_complete.cpp \
../example/cluster/example_cluster_linkage_single.cpp \
../example/cluster/example_cluster_linkage_total.cpp \
../example/cluster/example_cluster_node.cpp \
../example/cluster/example_cluster_node_colorer.cpp \
../example/cluster/example_cluster_node_description_average.cpp \
../example/cluster/example_cluster_node_description_from_file.cpp \
../example/cluster/example_cluster_output_centers.cpp \
../example/cluster/example_cluster_output_classes.cpp \
../example/cluster/example_cluster_output_pymol.cpp \
../example/cluster/example_cluster_output_pymol_label_small_molecule.cpp \
../example/cluster/example_cluster_output_pymol_label_string.cpp \
../example/cluster/example_cluster_output_rows.cpp 

OBJS += \
./example/cluster/example_cluster.o \
./example/cluster/example_cluster_dendrogram.o \
./example/cluster/example_cluster_distances_euclidean.o \
./example/cluster/example_cluster_distances_stored.o \
./example/cluster/example_cluster_input_classes.o \
./example/cluster/example_cluster_input_pairwise_list.o \
./example/cluster/example_cluster_input_table.o \
./example/cluster/example_cluster_linkage_average.o \
./example/cluster/example_cluster_linkage_classes.o \
./example/cluster/example_cluster_linkage_complete.o \
./example/cluster/example_cluster_linkage_single.o \
./example/cluster/example_cluster_linkage_total.o \
./example/cluster/example_cluster_node.o \
./example/cluster/example_cluster_node_colorer.o \
./example/cluster/example_cluster_node_description_average.o \
./example/cluster/example_cluster_node_description_from_file.o \
./example/cluster/example_cluster_output_centers.o \
./example/cluster/example_cluster_output_classes.o \
./example/cluster/example_cluster_output_pymol.o \
./example/cluster/example_cluster_output_pymol_label_small_molecule.o \
./example/cluster/example_cluster_output_pymol_label_string.o \
./example/cluster/example_cluster_output_rows.o 

CPP_DEPS += \
./example/cluster/example_cluster.d \
./example/cluster/example_cluster_dendrogram.d \
./example/cluster/example_cluster_distances_euclidean.d \
./example/cluster/example_cluster_distances_stored.d \
./example/cluster/example_cluster_input_classes.d \
./example/cluster/example_cluster_input_pairwise_list.d \
./example/cluster/example_cluster_input_table.d \
./example/cluster/example_cluster_linkage_average.d \
./example/cluster/example_cluster_linkage_classes.d \
./example/cluster/example_cluster_linkage_complete.d \
./example/cluster/example_cluster_linkage_single.d \
./example/cluster/example_cluster_linkage_total.d \
./example/cluster/example_cluster_node.d \
./example/cluster/example_cluster_node_colorer.d \
./example/cluster/example_cluster_node_description_average.d \
./example/cluster/example_cluster_node_description_from_file.d \
./example/cluster/example_cluster_output_centers.d \
./example/cluster/example_cluster_output_classes.d \
./example/cluster/example_cluster_output_pymol.d \
./example/cluster/example_cluster_output_pymol_label_small_molecule.d \
./example/cluster/example_cluster_output_pymol_label_string.d \
./example/cluster/example_cluster_output_rows.d 


# Each subdirectory must supply rules for building sources it contributes
example/cluster/%.o: ../example/cluster/%.cpp example/cluster/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



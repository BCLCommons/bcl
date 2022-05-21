################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/cluster/example_cluster.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_dendrogram.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_distances_euclidean.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_distances_stored.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_input_classes.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_input_pairwise_list.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_input_table.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_linkage_average.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_linkage_classes.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_linkage_complete.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_linkage_single.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_linkage_total.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_node.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_node_colorer.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_node_description_average.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_node_description_from_file.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_centers.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_classes.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_pymol.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_small_molecule.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_string.cpp \
../bcl_test_release/bcl/example/cluster/example_cluster_output_rows.cpp 

OBJS += \
./bcl_test_release/bcl/example/cluster/example_cluster.o \
./bcl_test_release/bcl/example/cluster/example_cluster_dendrogram.o \
./bcl_test_release/bcl/example/cluster/example_cluster_distances_euclidean.o \
./bcl_test_release/bcl/example/cluster/example_cluster_distances_stored.o \
./bcl_test_release/bcl/example/cluster/example_cluster_input_classes.o \
./bcl_test_release/bcl/example/cluster/example_cluster_input_pairwise_list.o \
./bcl_test_release/bcl/example/cluster/example_cluster_input_table.o \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_average.o \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_classes.o \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_complete.o \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_single.o \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_total.o \
./bcl_test_release/bcl/example/cluster/example_cluster_node.o \
./bcl_test_release/bcl/example/cluster/example_cluster_node_colorer.o \
./bcl_test_release/bcl/example/cluster/example_cluster_node_description_average.o \
./bcl_test_release/bcl/example/cluster/example_cluster_node_description_from_file.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_centers.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_classes.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_small_molecule.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_string.o \
./bcl_test_release/bcl/example/cluster/example_cluster_output_rows.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/cluster/example_cluster.d \
./bcl_test_release/bcl/example/cluster/example_cluster_dendrogram.d \
./bcl_test_release/bcl/example/cluster/example_cluster_distances_euclidean.d \
./bcl_test_release/bcl/example/cluster/example_cluster_distances_stored.d \
./bcl_test_release/bcl/example/cluster/example_cluster_input_classes.d \
./bcl_test_release/bcl/example/cluster/example_cluster_input_pairwise_list.d \
./bcl_test_release/bcl/example/cluster/example_cluster_input_table.d \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_average.d \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_classes.d \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_complete.d \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_single.d \
./bcl_test_release/bcl/example/cluster/example_cluster_linkage_total.d \
./bcl_test_release/bcl/example/cluster/example_cluster_node.d \
./bcl_test_release/bcl/example/cluster/example_cluster_node_colorer.d \
./bcl_test_release/bcl/example/cluster/example_cluster_node_description_average.d \
./bcl_test_release/bcl/example/cluster/example_cluster_node_description_from_file.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_centers.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_classes.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_small_molecule.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_pymol_label_string.d \
./bcl_test_release/bcl/example/cluster/example_cluster_output_rows.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/cluster/%.o: ../bcl_test_release/bcl/example/cluster/%.cpp bcl_test_release/bcl/example/cluster/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



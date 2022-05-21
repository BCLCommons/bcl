################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/graph/example_graph.cpp \
../bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism.cpp \
../bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism_base.cpp \
../bcl_test_release/bcl/example/graph/example_graph_const_graph.cpp \
../bcl_test_release/bcl/example/graph/example_graph_csi_substructure.cpp \
../bcl_test_release/bcl/example/graph/example_graph_edge_cover_ring_perception.cpp \
../bcl_test_release/bcl/example/graph/example_graph_edge_with_data.cpp \
../bcl_test_release/bcl/example/graph/example_graph_exhaustive_ring_perception.cpp \
../bcl_test_release/bcl/example/graph/example_graph_graph_with_data.cpp \
../bcl_test_release/bcl/example/graph/example_graph_path.cpp \
../bcl_test_release/bcl/example/graph/example_graph_ring.cpp \
../bcl_test_release/bcl/example/graph/example_graph_subgraph.cpp \
../bcl_test_release/bcl/example/graph/example_graph_subgraph_isomorphism.cpp \
../bcl_test_release/bcl/example/graph/example_graph_tree_node.cpp \
../bcl_test_release/bcl/example/graph/example_graph_undirected_edge.cpp \
../bcl_test_release/bcl/example/graph/example_graph_vertex_with_data.cpp 

OBJS += \
./bcl_test_release/bcl/example/graph/example_graph.o \
./bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism.o \
./bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism_base.o \
./bcl_test_release/bcl/example/graph/example_graph_const_graph.o \
./bcl_test_release/bcl/example/graph/example_graph_csi_substructure.o \
./bcl_test_release/bcl/example/graph/example_graph_edge_cover_ring_perception.o \
./bcl_test_release/bcl/example/graph/example_graph_edge_with_data.o \
./bcl_test_release/bcl/example/graph/example_graph_exhaustive_ring_perception.o \
./bcl_test_release/bcl/example/graph/example_graph_graph_with_data.o \
./bcl_test_release/bcl/example/graph/example_graph_path.o \
./bcl_test_release/bcl/example/graph/example_graph_ring.o \
./bcl_test_release/bcl/example/graph/example_graph_subgraph.o \
./bcl_test_release/bcl/example/graph/example_graph_subgraph_isomorphism.o \
./bcl_test_release/bcl/example/graph/example_graph_tree_node.o \
./bcl_test_release/bcl/example/graph/example_graph_undirected_edge.o \
./bcl_test_release/bcl/example/graph/example_graph_vertex_with_data.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/graph/example_graph.d \
./bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism.d \
./bcl_test_release/bcl/example/graph/example_graph_common_subgraph_isomorphism_base.d \
./bcl_test_release/bcl/example/graph/example_graph_const_graph.d \
./bcl_test_release/bcl/example/graph/example_graph_csi_substructure.d \
./bcl_test_release/bcl/example/graph/example_graph_edge_cover_ring_perception.d \
./bcl_test_release/bcl/example/graph/example_graph_edge_with_data.d \
./bcl_test_release/bcl/example/graph/example_graph_exhaustive_ring_perception.d \
./bcl_test_release/bcl/example/graph/example_graph_graph_with_data.d \
./bcl_test_release/bcl/example/graph/example_graph_path.d \
./bcl_test_release/bcl/example/graph/example_graph_ring.d \
./bcl_test_release/bcl/example/graph/example_graph_subgraph.d \
./bcl_test_release/bcl/example/graph/example_graph_subgraph_isomorphism.d \
./bcl_test_release/bcl/example/graph/example_graph_tree_node.d \
./bcl_test_release/bcl/example/graph/example_graph_undirected_edge.d \
./bcl_test_release/bcl/example/graph/example_graph_vertex_with_data.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/graph/%.o: ../bcl_test_release/bcl/example/graph/%.cpp bcl_test_release/bcl/example/graph/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



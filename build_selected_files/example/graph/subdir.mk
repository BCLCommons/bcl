################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/graph/example_graph.cpp \
../example/graph/example_graph_common_subgraph_isomorphism.cpp \
../example/graph/example_graph_common_subgraph_isomorphism_base.cpp \
../example/graph/example_graph_const_graph.cpp \
../example/graph/example_graph_csi_substructure.cpp \
../example/graph/example_graph_edge_cover_ring_perception.cpp \
../example/graph/example_graph_edge_with_data.cpp \
../example/graph/example_graph_exhaustive_ring_perception.cpp \
../example/graph/example_graph_graph_with_data.cpp \
../example/graph/example_graph_path.cpp \
../example/graph/example_graph_ring.cpp \
../example/graph/example_graph_subgraph.cpp \
../example/graph/example_graph_subgraph_isomorphism.cpp \
../example/graph/example_graph_tree_node.cpp \
../example/graph/example_graph_undirected_edge.cpp \
../example/graph/example_graph_vertex_with_data.cpp 

OBJS += \
./example/graph/example_graph.o \
./example/graph/example_graph_common_subgraph_isomorphism.o \
./example/graph/example_graph_common_subgraph_isomorphism_base.o \
./example/graph/example_graph_const_graph.o \
./example/graph/example_graph_csi_substructure.o \
./example/graph/example_graph_edge_cover_ring_perception.o \
./example/graph/example_graph_edge_with_data.o \
./example/graph/example_graph_exhaustive_ring_perception.o \
./example/graph/example_graph_graph_with_data.o \
./example/graph/example_graph_path.o \
./example/graph/example_graph_ring.o \
./example/graph/example_graph_subgraph.o \
./example/graph/example_graph_subgraph_isomorphism.o \
./example/graph/example_graph_tree_node.o \
./example/graph/example_graph_undirected_edge.o \
./example/graph/example_graph_vertex_with_data.o 

CPP_DEPS += \
./example/graph/example_graph.d \
./example/graph/example_graph_common_subgraph_isomorphism.d \
./example/graph/example_graph_common_subgraph_isomorphism_base.d \
./example/graph/example_graph_const_graph.d \
./example/graph/example_graph_csi_substructure.d \
./example/graph/example_graph_edge_cover_ring_perception.d \
./example/graph/example_graph_edge_with_data.d \
./example/graph/example_graph_exhaustive_ring_perception.d \
./example/graph/example_graph_graph_with_data.d \
./example/graph/example_graph_path.d \
./example/graph/example_graph_ring.d \
./example/graph/example_graph_subgraph.d \
./example/graph/example_graph_subgraph_isomorphism.d \
./example/graph/example_graph_tree_node.d \
./example/graph/example_graph_undirected_edge.d \
./example/graph/example_graph_vertex_with_data.d 


# Each subdirectory must supply rules for building sources it contributes
example/graph/%.o: ../example/graph/%.cpp example/graph/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



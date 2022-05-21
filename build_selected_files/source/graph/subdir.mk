################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/graph/bcl_graph.cpp \
../source/graph/bcl_graph_common_subgraph_isomorphism_base.cpp \
../source/graph/bcl_graph_edge_cover_ring_perception.cpp \
../source/graph/bcl_graph_exhaustive_ring_perception.cpp \
../source/graph/bcl_graph_path.cpp \
../source/graph/bcl_graph_ring.cpp \
../source/graph/bcl_graph_subgraph_isomorphism_base.cpp 

OBJS += \
./source/graph/bcl_graph.o \
./source/graph/bcl_graph_common_subgraph_isomorphism_base.o \
./source/graph/bcl_graph_edge_cover_ring_perception.o \
./source/graph/bcl_graph_exhaustive_ring_perception.o \
./source/graph/bcl_graph_path.o \
./source/graph/bcl_graph_ring.o \
./source/graph/bcl_graph_subgraph_isomorphism_base.o 

CPP_DEPS += \
./source/graph/bcl_graph.d \
./source/graph/bcl_graph_common_subgraph_isomorphism_base.d \
./source/graph/bcl_graph_edge_cover_ring_perception.d \
./source/graph/bcl_graph_exhaustive_ring_perception.d \
./source/graph/bcl_graph_path.d \
./source/graph/bcl_graph_ring.d \
./source/graph/bcl_graph_subgraph_isomorphism_base.d 


# Each subdirectory must supply rules for building sources it contributes
source/graph/%.o: ../source/graph/%.cpp source/graph/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



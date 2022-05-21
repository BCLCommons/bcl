################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/coord/example_coord.cpp \
../example/coord/example_coord_axes.cpp \
../example/coord/example_coord_cyclic_coordinate_descent.cpp \
../example/coord/example_coord_cylinder_coordinates.cpp \
../example/coord/example_coord_geometric_hash_storage_classes.cpp \
../example/coord/example_coord_geometric_hash_storage_hash_map.cpp \
../example/coord/example_coord_geometric_hashing.cpp \
../example/coord/example_coord_line_segment_2d.cpp \
../example/coord/example_coord_line_segment_3d.cpp \
../example/coord/example_coord_moment_of_inertia.cpp \
../example/coord/example_coord_movable_eccentric.cpp \
../example/coord/example_coord_move_combine.cpp \
../example/coord/example_coord_move_rotate_defined.cpp \
../example/coord/example_coord_move_rotate_random.cpp \
../example/coord/example_coord_move_rotate_random_external_reference.cpp \
../example/coord/example_coord_move_transform_random.cpp \
../example/coord/example_coord_move_translate_defined.cpp \
../example/coord/example_coord_move_translate_external_axis.cpp \
../example/coord/example_coord_move_translate_random.cpp \
../example/coord/example_coord_point_cloud.cpp \
../example/coord/example_coord_point_to_key_classes.cpp \
../example/coord/example_coord_point_to_key_spherical_radius.cpp \
../example/coord/example_coord_polygon.cpp \
../example/coord/example_coord_sphere.cpp 

OBJS += \
./example/coord/example_coord.o \
./example/coord/example_coord_axes.o \
./example/coord/example_coord_cyclic_coordinate_descent.o \
./example/coord/example_coord_cylinder_coordinates.o \
./example/coord/example_coord_geometric_hash_storage_classes.o \
./example/coord/example_coord_geometric_hash_storage_hash_map.o \
./example/coord/example_coord_geometric_hashing.o \
./example/coord/example_coord_line_segment_2d.o \
./example/coord/example_coord_line_segment_3d.o \
./example/coord/example_coord_moment_of_inertia.o \
./example/coord/example_coord_movable_eccentric.o \
./example/coord/example_coord_move_combine.o \
./example/coord/example_coord_move_rotate_defined.o \
./example/coord/example_coord_move_rotate_random.o \
./example/coord/example_coord_move_rotate_random_external_reference.o \
./example/coord/example_coord_move_transform_random.o \
./example/coord/example_coord_move_translate_defined.o \
./example/coord/example_coord_move_translate_external_axis.o \
./example/coord/example_coord_move_translate_random.o \
./example/coord/example_coord_point_cloud.o \
./example/coord/example_coord_point_to_key_classes.o \
./example/coord/example_coord_point_to_key_spherical_radius.o \
./example/coord/example_coord_polygon.o \
./example/coord/example_coord_sphere.o 

CPP_DEPS += \
./example/coord/example_coord.d \
./example/coord/example_coord_axes.d \
./example/coord/example_coord_cyclic_coordinate_descent.d \
./example/coord/example_coord_cylinder_coordinates.d \
./example/coord/example_coord_geometric_hash_storage_classes.d \
./example/coord/example_coord_geometric_hash_storage_hash_map.d \
./example/coord/example_coord_geometric_hashing.d \
./example/coord/example_coord_line_segment_2d.d \
./example/coord/example_coord_line_segment_3d.d \
./example/coord/example_coord_moment_of_inertia.d \
./example/coord/example_coord_movable_eccentric.d \
./example/coord/example_coord_move_combine.d \
./example/coord/example_coord_move_rotate_defined.d \
./example/coord/example_coord_move_rotate_random.d \
./example/coord/example_coord_move_rotate_random_external_reference.d \
./example/coord/example_coord_move_transform_random.d \
./example/coord/example_coord_move_translate_defined.d \
./example/coord/example_coord_move_translate_external_axis.d \
./example/coord/example_coord_move_translate_random.d \
./example/coord/example_coord_point_cloud.d \
./example/coord/example_coord_point_to_key_classes.d \
./example/coord/example_coord_point_to_key_spherical_radius.d \
./example/coord/example_coord_polygon.d \
./example/coord/example_coord_sphere.d 


# Each subdirectory must supply rules for building sources it contributes
example/coord/%.o: ../example/coord/%.cpp example/coord/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



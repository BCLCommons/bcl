################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/coord/bcl_coord.cpp \
../source/coord/bcl_coord_axes.cpp \
../source/coord/bcl_coord_cyclic_coordinate_descent.cpp \
../source/coord/bcl_coord_cylinder_coordinates.cpp \
../source/coord/bcl_coord_geometric_hash_storage_classes.cpp \
../source/coord/bcl_coord_geometric_hash_storage_hash_map.cpp \
../source/coord/bcl_coord_geometric_hash_storage_interface.cpp \
../source/coord/bcl_coord_geometric_hashing.cpp \
../source/coord/bcl_coord_geometry_interface.cpp \
../source/coord/bcl_coord_line_segment_2d.cpp \
../source/coord/bcl_coord_line_segment_3d.cpp \
../source/coord/bcl_coord_moment_of_inertia.cpp \
../source/coord/bcl_coord_movable_eccentric.cpp \
../source/coord/bcl_coord_movable_interface.cpp \
../source/coord/bcl_coord_move_combine.cpp \
../source/coord/bcl_coord_move_rotate_defined.cpp \
../source/coord/bcl_coord_move_rotate_random.cpp \
../source/coord/bcl_coord_move_rotate_random_external_reference.cpp \
../source/coord/bcl_coord_move_transform_random.cpp \
../source/coord/bcl_coord_move_translate_defined.cpp \
../source/coord/bcl_coord_move_translate_external_axis.cpp \
../source/coord/bcl_coord_move_translate_random.cpp \
../source/coord/bcl_coord_orientation_interface.cpp \
../source/coord/bcl_coord_point_cloud.cpp \
../source/coord/bcl_coord_point_to_key_cartesian.cpp \
../source/coord/bcl_coord_point_to_key_classes.cpp \
../source/coord/bcl_coord_point_to_key_interface.cpp \
../source/coord/bcl_coord_point_to_key_spherical.cpp \
../source/coord/bcl_coord_point_to_key_spherical_radius.cpp \
../source/coord/bcl_coord_polygon.cpp \
../source/coord/bcl_coord_sphere.cpp 

OBJS += \
./source/coord/bcl_coord.o \
./source/coord/bcl_coord_axes.o \
./source/coord/bcl_coord_cyclic_coordinate_descent.o \
./source/coord/bcl_coord_cylinder_coordinates.o \
./source/coord/bcl_coord_geometric_hash_storage_classes.o \
./source/coord/bcl_coord_geometric_hash_storage_hash_map.o \
./source/coord/bcl_coord_geometric_hash_storage_interface.o \
./source/coord/bcl_coord_geometric_hashing.o \
./source/coord/bcl_coord_geometry_interface.o \
./source/coord/bcl_coord_line_segment_2d.o \
./source/coord/bcl_coord_line_segment_3d.o \
./source/coord/bcl_coord_moment_of_inertia.o \
./source/coord/bcl_coord_movable_eccentric.o \
./source/coord/bcl_coord_movable_interface.o \
./source/coord/bcl_coord_move_combine.o \
./source/coord/bcl_coord_move_rotate_defined.o \
./source/coord/bcl_coord_move_rotate_random.o \
./source/coord/bcl_coord_move_rotate_random_external_reference.o \
./source/coord/bcl_coord_move_transform_random.o \
./source/coord/bcl_coord_move_translate_defined.o \
./source/coord/bcl_coord_move_translate_external_axis.o \
./source/coord/bcl_coord_move_translate_random.o \
./source/coord/bcl_coord_orientation_interface.o \
./source/coord/bcl_coord_point_cloud.o \
./source/coord/bcl_coord_point_to_key_cartesian.o \
./source/coord/bcl_coord_point_to_key_classes.o \
./source/coord/bcl_coord_point_to_key_interface.o \
./source/coord/bcl_coord_point_to_key_spherical.o \
./source/coord/bcl_coord_point_to_key_spherical_radius.o \
./source/coord/bcl_coord_polygon.o \
./source/coord/bcl_coord_sphere.o 

CPP_DEPS += \
./source/coord/bcl_coord.d \
./source/coord/bcl_coord_axes.d \
./source/coord/bcl_coord_cyclic_coordinate_descent.d \
./source/coord/bcl_coord_cylinder_coordinates.d \
./source/coord/bcl_coord_geometric_hash_storage_classes.d \
./source/coord/bcl_coord_geometric_hash_storage_hash_map.d \
./source/coord/bcl_coord_geometric_hash_storage_interface.d \
./source/coord/bcl_coord_geometric_hashing.d \
./source/coord/bcl_coord_geometry_interface.d \
./source/coord/bcl_coord_line_segment_2d.d \
./source/coord/bcl_coord_line_segment_3d.d \
./source/coord/bcl_coord_moment_of_inertia.d \
./source/coord/bcl_coord_movable_eccentric.d \
./source/coord/bcl_coord_movable_interface.d \
./source/coord/bcl_coord_move_combine.d \
./source/coord/bcl_coord_move_rotate_defined.d \
./source/coord/bcl_coord_move_rotate_random.d \
./source/coord/bcl_coord_move_rotate_random_external_reference.d \
./source/coord/bcl_coord_move_transform_random.d \
./source/coord/bcl_coord_move_translate_defined.d \
./source/coord/bcl_coord_move_translate_external_axis.d \
./source/coord/bcl_coord_move_translate_random.d \
./source/coord/bcl_coord_orientation_interface.d \
./source/coord/bcl_coord_point_cloud.d \
./source/coord/bcl_coord_point_to_key_cartesian.d \
./source/coord/bcl_coord_point_to_key_classes.d \
./source/coord/bcl_coord_point_to_key_interface.d \
./source/coord/bcl_coord_point_to_key_spherical.d \
./source/coord/bcl_coord_point_to_key_spherical_radius.d \
./source/coord/bcl_coord_polygon.d \
./source/coord/bcl_coord_sphere.d 


# Each subdirectory must supply rules for building sources it contributes
source/coord/%.o: ../source/coord/%.cpp source/coord/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



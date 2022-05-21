################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/apps/molecule/example_app_align_to_scaffold.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_build_scaffold_library.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_compare.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_coordinates.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_filter.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_properties.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_reorder.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_split.cpp \
../bcl_test_release/bcl/example/apps/molecule/example_app_molecule_unique.cpp 

OBJS += \
./bcl_test_release/bcl/example/apps/molecule/example_app_align_to_scaffold.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_build_scaffold_library.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_compare.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_coordinates.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_filter.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_properties.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_reorder.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_split.o \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_unique.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/apps/molecule/example_app_align_to_scaffold.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_build_scaffold_library.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_compare.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_coordinates.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_filter.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_properties.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_reorder.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_split.d \
./bcl_test_release/bcl/example/apps/molecule/example_app_molecule_unique.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/apps/molecule/%.o: ../bcl_test_release/bcl/example/apps/molecule/%.cpp bcl_test_release/bcl/example/apps/molecule/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/apps/molecule/example_app_align_to_scaffold.cpp \
../example/apps/molecule/example_app_build_scaffold_library.cpp \
../example/apps/molecule/example_app_molecule_compare.cpp \
../example/apps/molecule/example_app_molecule_coordinates.cpp \
../example/apps/molecule/example_app_molecule_filter.cpp \
../example/apps/molecule/example_app_molecule_properties.cpp \
../example/apps/molecule/example_app_molecule_reorder.cpp \
../example/apps/molecule/example_app_molecule_split.cpp \
../example/apps/molecule/example_app_molecule_unique.cpp 

OBJS += \
./example/apps/molecule/example_app_align_to_scaffold.o \
./example/apps/molecule/example_app_build_scaffold_library.o \
./example/apps/molecule/example_app_molecule_compare.o \
./example/apps/molecule/example_app_molecule_coordinates.o \
./example/apps/molecule/example_app_molecule_filter.o \
./example/apps/molecule/example_app_molecule_properties.o \
./example/apps/molecule/example_app_molecule_reorder.o \
./example/apps/molecule/example_app_molecule_split.o \
./example/apps/molecule/example_app_molecule_unique.o 

CPP_DEPS += \
./example/apps/molecule/example_app_align_to_scaffold.d \
./example/apps/molecule/example_app_build_scaffold_library.d \
./example/apps/molecule/example_app_molecule_compare.d \
./example/apps/molecule/example_app_molecule_coordinates.d \
./example/apps/molecule/example_app_molecule_filter.d \
./example/apps/molecule/example_app_molecule_properties.d \
./example/apps/molecule/example_app_molecule_reorder.d \
./example/apps/molecule/example_app_molecule_split.d \
./example/apps/molecule/example_app_molecule_unique.d 


# Each subdirectory must supply rules for building sources it contributes
example/apps/molecule/%.o: ../example/apps/molecule/%.cpp example/apps/molecule/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



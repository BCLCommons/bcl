################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/sdf/example_sdf.cpp \
../example/sdf/example_sdf_atom_info.cpp \
../example/sdf/example_sdf_bond_info.cpp \
../example/sdf/example_sdf_ctab_handler.cpp \
../example/sdf/example_sdf_factory.cpp \
../example/sdf/example_sdf_fragment_factory.cpp \
../example/sdf/example_sdf_mdl_entry_type_data.cpp \
../example/sdf/example_sdf_mdl_entry_types.cpp \
../example/sdf/example_sdf_mdl_handler.cpp \
../example/sdf/example_sdf_mdl_header.cpp \
../example/sdf/example_sdf_mdl_property.cpp \
../example/sdf/example_sdf_molfile_handler.cpp \
../example/sdf/example_sdf_rxn_factory.cpp \
../example/sdf/example_sdf_rxn_handler.cpp 

OBJS += \
./example/sdf/example_sdf.o \
./example/sdf/example_sdf_atom_info.o \
./example/sdf/example_sdf_bond_info.o \
./example/sdf/example_sdf_ctab_handler.o \
./example/sdf/example_sdf_factory.o \
./example/sdf/example_sdf_fragment_factory.o \
./example/sdf/example_sdf_mdl_entry_type_data.o \
./example/sdf/example_sdf_mdl_entry_types.o \
./example/sdf/example_sdf_mdl_handler.o \
./example/sdf/example_sdf_mdl_header.o \
./example/sdf/example_sdf_mdl_property.o \
./example/sdf/example_sdf_molfile_handler.o \
./example/sdf/example_sdf_rxn_factory.o \
./example/sdf/example_sdf_rxn_handler.o 

CPP_DEPS += \
./example/sdf/example_sdf.d \
./example/sdf/example_sdf_atom_info.d \
./example/sdf/example_sdf_bond_info.d \
./example/sdf/example_sdf_ctab_handler.d \
./example/sdf/example_sdf_factory.d \
./example/sdf/example_sdf_fragment_factory.d \
./example/sdf/example_sdf_mdl_entry_type_data.d \
./example/sdf/example_sdf_mdl_entry_types.d \
./example/sdf/example_sdf_mdl_handler.d \
./example/sdf/example_sdf_mdl_header.d \
./example/sdf/example_sdf_mdl_property.d \
./example/sdf/example_sdf_molfile_handler.d \
./example/sdf/example_sdf_rxn_factory.d \
./example/sdf/example_sdf_rxn_handler.d 


# Each subdirectory must supply rules for building sources it contributes
example/sdf/%.o: ../example/sdf/%.cpp example/sdf/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



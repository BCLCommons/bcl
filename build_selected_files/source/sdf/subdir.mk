################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/sdf/bcl_sdf.cpp \
../source/sdf/bcl_sdf_atom_info.cpp \
../source/sdf/bcl_sdf_bond_info.cpp \
../source/sdf/bcl_sdf_ctab_handler.cpp \
../source/sdf/bcl_sdf_factory.cpp \
../source/sdf/bcl_sdf_fragment_factory.cpp \
../source/sdf/bcl_sdf_mdl_entry_type_data.cpp \
../source/sdf/bcl_sdf_mdl_entry_types.cpp \
../source/sdf/bcl_sdf_mdl_handler.cpp \
../source/sdf/bcl_sdf_mdl_header.cpp \
../source/sdf/bcl_sdf_mdl_line_types.cpp \
../source/sdf/bcl_sdf_mdl_property.cpp \
../source/sdf/bcl_sdf_molecule_reading_pref.cpp \
../source/sdf/bcl_sdf_molfile_handler.cpp \
../source/sdf/bcl_sdf_rxn_factory.cpp \
../source/sdf/bcl_sdf_rxn_handler.cpp 

OBJS += \
./source/sdf/bcl_sdf.o \
./source/sdf/bcl_sdf_atom_info.o \
./source/sdf/bcl_sdf_bond_info.o \
./source/sdf/bcl_sdf_ctab_handler.o \
./source/sdf/bcl_sdf_factory.o \
./source/sdf/bcl_sdf_fragment_factory.o \
./source/sdf/bcl_sdf_mdl_entry_type_data.o \
./source/sdf/bcl_sdf_mdl_entry_types.o \
./source/sdf/bcl_sdf_mdl_handler.o \
./source/sdf/bcl_sdf_mdl_header.o \
./source/sdf/bcl_sdf_mdl_line_types.o \
./source/sdf/bcl_sdf_mdl_property.o \
./source/sdf/bcl_sdf_molecule_reading_pref.o \
./source/sdf/bcl_sdf_molfile_handler.o \
./source/sdf/bcl_sdf_rxn_factory.o \
./source/sdf/bcl_sdf_rxn_handler.o 

CPP_DEPS += \
./source/sdf/bcl_sdf.d \
./source/sdf/bcl_sdf_atom_info.d \
./source/sdf/bcl_sdf_bond_info.d \
./source/sdf/bcl_sdf_ctab_handler.d \
./source/sdf/bcl_sdf_factory.d \
./source/sdf/bcl_sdf_fragment_factory.d \
./source/sdf/bcl_sdf_mdl_entry_type_data.d \
./source/sdf/bcl_sdf_mdl_entry_types.d \
./source/sdf/bcl_sdf_mdl_handler.d \
./source/sdf/bcl_sdf_mdl_header.d \
./source/sdf/bcl_sdf_mdl_line_types.d \
./source/sdf/bcl_sdf_mdl_property.d \
./source/sdf/bcl_sdf_molecule_reading_pref.d \
./source/sdf/bcl_sdf_molfile_handler.d \
./source/sdf/bcl_sdf_rxn_factory.d \
./source/sdf/bcl_sdf_rxn_handler.d 


# Each subdirectory must supply rules for building sources it contributes
source/sdf/%.o: ../source/sdf/%.cpp source/sdf/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



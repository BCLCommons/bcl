################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/sdf/example_sdf.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_atom_info.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_bond_info.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_ctab_handler.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_factory.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_fragment_factory.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_type_data.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_types.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_mdl_handler.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_mdl_header.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_mdl_property.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_molfile_handler.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_rxn_factory.cpp \
../bcl_test_release/bcl/example/sdf/example_sdf_rxn_handler.cpp 

OBJS += \
./bcl_test_release/bcl/example/sdf/example_sdf.o \
./bcl_test_release/bcl/example/sdf/example_sdf_atom_info.o \
./bcl_test_release/bcl/example/sdf/example_sdf_bond_info.o \
./bcl_test_release/bcl/example/sdf/example_sdf_ctab_handler.o \
./bcl_test_release/bcl/example/sdf/example_sdf_factory.o \
./bcl_test_release/bcl/example/sdf/example_sdf_fragment_factory.o \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_type_data.o \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_types.o \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_handler.o \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_header.o \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_property.o \
./bcl_test_release/bcl/example/sdf/example_sdf_molfile_handler.o \
./bcl_test_release/bcl/example/sdf/example_sdf_rxn_factory.o \
./bcl_test_release/bcl/example/sdf/example_sdf_rxn_handler.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/sdf/example_sdf.d \
./bcl_test_release/bcl/example/sdf/example_sdf_atom_info.d \
./bcl_test_release/bcl/example/sdf/example_sdf_bond_info.d \
./bcl_test_release/bcl/example/sdf/example_sdf_ctab_handler.d \
./bcl_test_release/bcl/example/sdf/example_sdf_factory.d \
./bcl_test_release/bcl/example/sdf/example_sdf_fragment_factory.d \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_type_data.d \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_entry_types.d \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_handler.d \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_header.d \
./bcl_test_release/bcl/example/sdf/example_sdf_mdl_property.d \
./bcl_test_release/bcl/example/sdf/example_sdf_molfile_handler.d \
./bcl_test_release/bcl/example/sdf/example_sdf_rxn_factory.d \
./bcl_test_release/bcl/example/sdf/example_sdf_rxn_handler.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/sdf/%.o: ../bcl_test_release/bcl/example/sdf/%.cpp bcl_test_release/bcl/example/sdf/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



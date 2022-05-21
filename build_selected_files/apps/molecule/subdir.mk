################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/molecule/bcl_app_align_to_scaffold.cpp \
../apps/molecule/bcl_app_evolutionary_generator.cpp \
../apps/molecule/bcl_app_extract_mcs.cpp \
../apps/molecule/bcl_app_focused_library_design.cpp \
../apps/molecule/bcl_app_generate_rosetta_ncaa_instructions.cpp \
../apps/molecule/bcl_app_make_fld_scaffold.cpp \
../apps/molecule/bcl_app_make_grow_fragments.cpp \
../apps/molecule/bcl_app_map_params.cpp \
../apps/molecule/bcl_app_match_reactions.cpp \
../apps/molecule/bcl_app_molecule_atom_pairwise_potential.cpp \
../apps/molecule/bcl_app_molecule_compare.cpp \
../apps/molecule/bcl_app_molecule_coordinates.cpp \
../apps/molecule/bcl_app_molecule_extract_protein_pocket.cpp \
../apps/molecule/bcl_app_molecule_features.cpp \
../apps/molecule/bcl_app_molecule_filter.cpp \
../apps/molecule/bcl_app_molecule_fit.cpp \
../apps/molecule/bcl_app_molecule_multialign.cpp \
../apps/molecule/bcl_app_molecule_properties.cpp \
../apps/molecule/bcl_app_molecule_reorder.cpp \
../apps/molecule/bcl_app_molecule_split.cpp \
../apps/molecule/bcl_app_molecule_unique.cpp \
../apps/molecule/bcl_app_pharm_map.cpp \
../apps/molecule/bcl_app_set_sample_by_parts_atoms.cpp 

OBJS += \
./apps/molecule/bcl_app_align_to_scaffold.o \
./apps/molecule/bcl_app_evolutionary_generator.o \
./apps/molecule/bcl_app_extract_mcs.o \
./apps/molecule/bcl_app_focused_library_design.o \
./apps/molecule/bcl_app_generate_rosetta_ncaa_instructions.o \
./apps/molecule/bcl_app_make_fld_scaffold.o \
./apps/molecule/bcl_app_make_grow_fragments.o \
./apps/molecule/bcl_app_map_params.o \
./apps/molecule/bcl_app_match_reactions.o \
./apps/molecule/bcl_app_molecule_atom_pairwise_potential.o \
./apps/molecule/bcl_app_molecule_compare.o \
./apps/molecule/bcl_app_molecule_coordinates.o \
./apps/molecule/bcl_app_molecule_extract_protein_pocket.o \
./apps/molecule/bcl_app_molecule_features.o \
./apps/molecule/bcl_app_molecule_filter.o \
./apps/molecule/bcl_app_molecule_fit.o \
./apps/molecule/bcl_app_molecule_multialign.o \
./apps/molecule/bcl_app_molecule_properties.o \
./apps/molecule/bcl_app_molecule_reorder.o \
./apps/molecule/bcl_app_molecule_split.o \
./apps/molecule/bcl_app_molecule_unique.o \
./apps/molecule/bcl_app_pharm_map.o \
./apps/molecule/bcl_app_set_sample_by_parts_atoms.o 

CPP_DEPS += \
./apps/molecule/bcl_app_align_to_scaffold.d \
./apps/molecule/bcl_app_evolutionary_generator.d \
./apps/molecule/bcl_app_extract_mcs.d \
./apps/molecule/bcl_app_focused_library_design.d \
./apps/molecule/bcl_app_generate_rosetta_ncaa_instructions.d \
./apps/molecule/bcl_app_make_fld_scaffold.d \
./apps/molecule/bcl_app_make_grow_fragments.d \
./apps/molecule/bcl_app_map_params.d \
./apps/molecule/bcl_app_match_reactions.d \
./apps/molecule/bcl_app_molecule_atom_pairwise_potential.d \
./apps/molecule/bcl_app_molecule_compare.d \
./apps/molecule/bcl_app_molecule_coordinates.d \
./apps/molecule/bcl_app_molecule_extract_protein_pocket.d \
./apps/molecule/bcl_app_molecule_features.d \
./apps/molecule/bcl_app_molecule_filter.d \
./apps/molecule/bcl_app_molecule_fit.d \
./apps/molecule/bcl_app_molecule_multialign.d \
./apps/molecule/bcl_app_molecule_properties.d \
./apps/molecule/bcl_app_molecule_reorder.d \
./apps/molecule/bcl_app_molecule_split.d \
./apps/molecule/bcl_app_molecule_unique.d \
./apps/molecule/bcl_app_pharm_map.d \
./apps/molecule/bcl_app_set_sample_by_parts_atoms.d 


# Each subdirectory must supply rules for building sources it contributes
apps/molecule/%.o: ../apps/molecule/%.cpp apps/molecule/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



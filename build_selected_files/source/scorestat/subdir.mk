################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/scorestat/bcl_scorestat.cpp \
../source/scorestat/bcl_scorestat_aa_count.cpp \
../source/scorestat/bcl_scorestat_aa_distance.cpp \
../source/scorestat/bcl_scorestat_aa_distance_angle_contacts.cpp \
../source/scorestat/bcl_scorestat_aa_distance_matrix.cpp \
../source/scorestat/bcl_scorestat_contact_order.cpp \
../source/scorestat/bcl_scorestat_fold_template.cpp \
../source/scorestat/bcl_scorestat_loop_angle.cpp \
../source/scorestat/bcl_scorestat_loop_closure.cpp \
../source/scorestat/bcl_scorestat_loop_distance.cpp \
../source/scorestat/bcl_scorestat_neighbor_count.cpp \
../source/scorestat/bcl_scorestat_neighbor_vector.cpp \
../source/scorestat/bcl_scorestat_ols.cpp \
../source/scorestat/bcl_scorestat_phipsi.cpp \
../source/scorestat/bcl_scorestat_protein_model_packing.cpp \
../source/scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.cpp \
../source/scorestat/bcl_scorestat_radius_of_gyration.cpp \
../source/scorestat/bcl_scorestat_sheet_template.cpp \
../source/scorestat/bcl_scorestat_side_chain_distance.cpp \
../source/scorestat/bcl_scorestat_sse_count.cpp \
../source/scorestat/bcl_scorestat_sse_membrane_alignment.cpp \
../source/scorestat/bcl_scorestat_sse_packing.cpp \
../source/scorestat/bcl_scorestat_sspred_agreement.cpp \
../source/scorestat/bcl_scorestat_strand_alignment.cpp 

OBJS += \
./source/scorestat/bcl_scorestat.o \
./source/scorestat/bcl_scorestat_aa_count.o \
./source/scorestat/bcl_scorestat_aa_distance.o \
./source/scorestat/bcl_scorestat_aa_distance_angle_contacts.o \
./source/scorestat/bcl_scorestat_aa_distance_matrix.o \
./source/scorestat/bcl_scorestat_contact_order.o \
./source/scorestat/bcl_scorestat_fold_template.o \
./source/scorestat/bcl_scorestat_loop_angle.o \
./source/scorestat/bcl_scorestat_loop_closure.o \
./source/scorestat/bcl_scorestat_loop_distance.o \
./source/scorestat/bcl_scorestat_neighbor_count.o \
./source/scorestat/bcl_scorestat_neighbor_vector.o \
./source/scorestat/bcl_scorestat_ols.o \
./source/scorestat/bcl_scorestat_phipsi.o \
./source/scorestat/bcl_scorestat_protein_model_packing.o \
./source/scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.o \
./source/scorestat/bcl_scorestat_radius_of_gyration.o \
./source/scorestat/bcl_scorestat_sheet_template.o \
./source/scorestat/bcl_scorestat_side_chain_distance.o \
./source/scorestat/bcl_scorestat_sse_count.o \
./source/scorestat/bcl_scorestat_sse_membrane_alignment.o \
./source/scorestat/bcl_scorestat_sse_packing.o \
./source/scorestat/bcl_scorestat_sspred_agreement.o \
./source/scorestat/bcl_scorestat_strand_alignment.o 

CPP_DEPS += \
./source/scorestat/bcl_scorestat.d \
./source/scorestat/bcl_scorestat_aa_count.d \
./source/scorestat/bcl_scorestat_aa_distance.d \
./source/scorestat/bcl_scorestat_aa_distance_angle_contacts.d \
./source/scorestat/bcl_scorestat_aa_distance_matrix.d \
./source/scorestat/bcl_scorestat_contact_order.d \
./source/scorestat/bcl_scorestat_fold_template.d \
./source/scorestat/bcl_scorestat_loop_angle.d \
./source/scorestat/bcl_scorestat_loop_closure.d \
./source/scorestat/bcl_scorestat_loop_distance.d \
./source/scorestat/bcl_scorestat_neighbor_count.d \
./source/scorestat/bcl_scorestat_neighbor_vector.d \
./source/scorestat/bcl_scorestat_ols.d \
./source/scorestat/bcl_scorestat_phipsi.d \
./source/scorestat/bcl_scorestat_protein_model_packing.d \
./source/scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.d \
./source/scorestat/bcl_scorestat_radius_of_gyration.d \
./source/scorestat/bcl_scorestat_sheet_template.d \
./source/scorestat/bcl_scorestat_side_chain_distance.d \
./source/scorestat/bcl_scorestat_sse_count.d \
./source/scorestat/bcl_scorestat_sse_membrane_alignment.d \
./source/scorestat/bcl_scorestat_sse_packing.d \
./source/scorestat/bcl_scorestat_sspred_agreement.d \
./source/scorestat/bcl_scorestat_strand_alignment.d 


# Each subdirectory must supply rules for building sources it contributes
source/scorestat/%.o: ../source/scorestat/%.cpp source/scorestat/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/descriptor/bcl_app_descriptor_analyze.cpp \
../apps/descriptor/bcl_app_descriptor_convert_code_object_file.cpp \
../apps/descriptor/bcl_app_descriptor_dataset_similarity_measures.cpp \
../apps/descriptor/bcl_app_descriptor_generate_dataset.cpp \
../apps/descriptor/bcl_app_descriptor_generate_pca_eigenvectors.cpp \
../apps/descriptor/bcl_app_descriptor_refine_by_score.cpp \
../apps/descriptor/bcl_app_descriptor_score_dataset.cpp \
../apps/descriptor/bcl_app_descriptor_sequential_feature_selection.cpp 

OBJS += \
./apps/descriptor/bcl_app_descriptor_analyze.o \
./apps/descriptor/bcl_app_descriptor_convert_code_object_file.o \
./apps/descriptor/bcl_app_descriptor_dataset_similarity_measures.o \
./apps/descriptor/bcl_app_descriptor_generate_dataset.o \
./apps/descriptor/bcl_app_descriptor_generate_pca_eigenvectors.o \
./apps/descriptor/bcl_app_descriptor_refine_by_score.o \
./apps/descriptor/bcl_app_descriptor_score_dataset.o \
./apps/descriptor/bcl_app_descriptor_sequential_feature_selection.o 

CPP_DEPS += \
./apps/descriptor/bcl_app_descriptor_analyze.d \
./apps/descriptor/bcl_app_descriptor_convert_code_object_file.d \
./apps/descriptor/bcl_app_descriptor_dataset_similarity_measures.d \
./apps/descriptor/bcl_app_descriptor_generate_dataset.d \
./apps/descriptor/bcl_app_descriptor_generate_pca_eigenvectors.d \
./apps/descriptor/bcl_app_descriptor_refine_by_score.d \
./apps/descriptor/bcl_app_descriptor_score_dataset.d \
./apps/descriptor/bcl_app_descriptor_sequential_feature_selection.d 


# Each subdirectory must supply rules for building sources it contributes
apps/descriptor/%.o: ../apps/descriptor/%.cpp apps/descriptor/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



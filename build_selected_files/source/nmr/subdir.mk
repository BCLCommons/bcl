################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/nmr/bcl_nmr.cpp \
../source/nmr/bcl_nmr_rdc_container.cpp \
../source/nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.cpp \
../source/nmr/bcl_nmr_rosetta_noe_handler.cpp \
../source/nmr/bcl_nmr_rosetta_rdc_handler.cpp \
../source/nmr/bcl_nmr_signal.cpp \
../source/nmr/bcl_nmr_signal_1d.cpp \
../source/nmr/bcl_nmr_spectrum.cpp \
../source/nmr/bcl_nmr_star_noe_handler.cpp \
../source/nmr/bcl_nmr_star_rdc_handler.cpp \
../source/nmr/bcl_nmr_star_tag_categories.cpp \
../source/nmr/bcl_nmr_star_tag_category_data.cpp \
../source/nmr/bcl_nmr_star_tag_data.cpp \
../source/nmr/bcl_nmr_star_tags.cpp 

OBJS += \
./source/nmr/bcl_nmr.o \
./source/nmr/bcl_nmr_rdc_container.o \
./source/nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.o \
./source/nmr/bcl_nmr_rosetta_noe_handler.o \
./source/nmr/bcl_nmr_rosetta_rdc_handler.o \
./source/nmr/bcl_nmr_signal.o \
./source/nmr/bcl_nmr_signal_1d.o \
./source/nmr/bcl_nmr_spectrum.o \
./source/nmr/bcl_nmr_star_noe_handler.o \
./source/nmr/bcl_nmr_star_rdc_handler.o \
./source/nmr/bcl_nmr_star_tag_categories.o \
./source/nmr/bcl_nmr_star_tag_category_data.o \
./source/nmr/bcl_nmr_star_tag_data.o \
./source/nmr/bcl_nmr_star_tags.o 

CPP_DEPS += \
./source/nmr/bcl_nmr.d \
./source/nmr/bcl_nmr_rdc_container.d \
./source/nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.d \
./source/nmr/bcl_nmr_rosetta_noe_handler.d \
./source/nmr/bcl_nmr_rosetta_rdc_handler.d \
./source/nmr/bcl_nmr_signal.d \
./source/nmr/bcl_nmr_signal_1d.d \
./source/nmr/bcl_nmr_spectrum.d \
./source/nmr/bcl_nmr_star_noe_handler.d \
./source/nmr/bcl_nmr_star_rdc_handler.d \
./source/nmr/bcl_nmr_star_tag_categories.d \
./source/nmr/bcl_nmr_star_tag_category_data.d \
./source/nmr/bcl_nmr_star_tag_data.d \
./source/nmr/bcl_nmr_star_tags.d 


# Each subdirectory must supply rules for building sources it contributes
source/nmr/%.o: ../source/nmr/%.cpp source/nmr/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



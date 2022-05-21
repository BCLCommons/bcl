################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/nmr/example_nmr.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_rdc_container.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_residual_dipolar_coupling_least_square_deviation.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_rosetta_noe_handler.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_rosetta_rdc_handler.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_signal.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_signal_1d.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_spectrum.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_noe_handler.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_rdc_handler.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_tag_categories.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_tag_category_data.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_tag_data.cpp \
../bcl_test_release/bcl/example/nmr/example_nmr_star_tags.cpp 

OBJS += \
./bcl_test_release/bcl/example/nmr/example_nmr.o \
./bcl_test_release/bcl/example/nmr/example_nmr_rdc_container.o \
./bcl_test_release/bcl/example/nmr/example_nmr_residual_dipolar_coupling_least_square_deviation.o \
./bcl_test_release/bcl/example/nmr/example_nmr_rosetta_noe_handler.o \
./bcl_test_release/bcl/example/nmr/example_nmr_rosetta_rdc_handler.o \
./bcl_test_release/bcl/example/nmr/example_nmr_signal.o \
./bcl_test_release/bcl/example/nmr/example_nmr_signal_1d.o \
./bcl_test_release/bcl/example/nmr/example_nmr_spectrum.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_noe_handler.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_rdc_handler.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_categories.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_category_data.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_data.o \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tags.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/nmr/example_nmr.d \
./bcl_test_release/bcl/example/nmr/example_nmr_rdc_container.d \
./bcl_test_release/bcl/example/nmr/example_nmr_residual_dipolar_coupling_least_square_deviation.d \
./bcl_test_release/bcl/example/nmr/example_nmr_rosetta_noe_handler.d \
./bcl_test_release/bcl/example/nmr/example_nmr_rosetta_rdc_handler.d \
./bcl_test_release/bcl/example/nmr/example_nmr_signal.d \
./bcl_test_release/bcl/example/nmr/example_nmr_signal_1d.d \
./bcl_test_release/bcl/example/nmr/example_nmr_spectrum.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_noe_handler.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_rdc_handler.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_categories.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_category_data.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tag_data.d \
./bcl_test_release/bcl/example/nmr/example_nmr_star_tags.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/nmr/%.o: ../bcl_test_release/bcl/example/nmr/%.cpp bcl_test_release/bcl/example/nmr/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



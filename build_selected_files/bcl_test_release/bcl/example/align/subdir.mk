################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../bcl_test_release/bcl/example/align/example_align.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_dp.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_dynamic_programming.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_merge.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_progressive.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_shift.cpp \
../bcl_test_release/bcl/example/align/example_align_aligner_wordbased.cpp \
../bcl_test_release/bcl/example/align/example_align_alignment_hit.cpp \
../bcl_test_release/bcl/example/align/example_align_alignment_leaf.cpp \
../bcl_test_release/bcl/example/align/example_align_alignment_node.cpp \
../bcl_test_release/bcl/example/align/example_align_alignment_word.cpp \
../bcl_test_release/bcl/example/align/example_align_assignment.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_blc.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_blocked.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_classes.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_fasta.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_pir.cpp \
../bcl_test_release/bcl/example/align/example_align_handler_standard.cpp \
../bcl_test_release/bcl/example/align/example_align_multiple_aligner_classes.cpp \
../bcl_test_release/bcl/example/align/example_align_pairwise_aligner_classes.cpp \
../bcl_test_release/bcl/example/align/example_align_sequence.cpp \
../bcl_test_release/bcl/example/align/example_align_word_generator_high_scoring.cpp \
../bcl_test_release/bcl/example/align/example_align_word_generator_subsequences.cpp 

OBJS += \
./bcl_test_release/bcl/example/align/example_align.o \
./bcl_test_release/bcl/example/align/example_align_aligner_dp.o \
./bcl_test_release/bcl/example/align/example_align_aligner_dynamic_programming.o \
./bcl_test_release/bcl/example/align/example_align_aligner_merge.o \
./bcl_test_release/bcl/example/align/example_align_aligner_progressive.o \
./bcl_test_release/bcl/example/align/example_align_aligner_shift.o \
./bcl_test_release/bcl/example/align/example_align_aligner_wordbased.o \
./bcl_test_release/bcl/example/align/example_align_alignment_hit.o \
./bcl_test_release/bcl/example/align/example_align_alignment_leaf.o \
./bcl_test_release/bcl/example/align/example_align_alignment_node.o \
./bcl_test_release/bcl/example/align/example_align_alignment_word.o \
./bcl_test_release/bcl/example/align/example_align_assignment.o \
./bcl_test_release/bcl/example/align/example_align_handler_blc.o \
./bcl_test_release/bcl/example/align/example_align_handler_blocked.o \
./bcl_test_release/bcl/example/align/example_align_handler_classes.o \
./bcl_test_release/bcl/example/align/example_align_handler_fasta.o \
./bcl_test_release/bcl/example/align/example_align_handler_pir.o \
./bcl_test_release/bcl/example/align/example_align_handler_standard.o \
./bcl_test_release/bcl/example/align/example_align_multiple_aligner_classes.o \
./bcl_test_release/bcl/example/align/example_align_pairwise_aligner_classes.o \
./bcl_test_release/bcl/example/align/example_align_sequence.o \
./bcl_test_release/bcl/example/align/example_align_word_generator_high_scoring.o \
./bcl_test_release/bcl/example/align/example_align_word_generator_subsequences.o 

CPP_DEPS += \
./bcl_test_release/bcl/example/align/example_align.d \
./bcl_test_release/bcl/example/align/example_align_aligner_dp.d \
./bcl_test_release/bcl/example/align/example_align_aligner_dynamic_programming.d \
./bcl_test_release/bcl/example/align/example_align_aligner_merge.d \
./bcl_test_release/bcl/example/align/example_align_aligner_progressive.d \
./bcl_test_release/bcl/example/align/example_align_aligner_shift.d \
./bcl_test_release/bcl/example/align/example_align_aligner_wordbased.d \
./bcl_test_release/bcl/example/align/example_align_alignment_hit.d \
./bcl_test_release/bcl/example/align/example_align_alignment_leaf.d \
./bcl_test_release/bcl/example/align/example_align_alignment_node.d \
./bcl_test_release/bcl/example/align/example_align_alignment_word.d \
./bcl_test_release/bcl/example/align/example_align_assignment.d \
./bcl_test_release/bcl/example/align/example_align_handler_blc.d \
./bcl_test_release/bcl/example/align/example_align_handler_blocked.d \
./bcl_test_release/bcl/example/align/example_align_handler_classes.d \
./bcl_test_release/bcl/example/align/example_align_handler_fasta.d \
./bcl_test_release/bcl/example/align/example_align_handler_pir.d \
./bcl_test_release/bcl/example/align/example_align_handler_standard.d \
./bcl_test_release/bcl/example/align/example_align_multiple_aligner_classes.d \
./bcl_test_release/bcl/example/align/example_align_pairwise_aligner_classes.d \
./bcl_test_release/bcl/example/align/example_align_sequence.d \
./bcl_test_release/bcl/example/align/example_align_word_generator_high_scoring.d \
./bcl_test_release/bcl/example/align/example_align_word_generator_subsequences.d 


# Each subdirectory must supply rules for building sources it contributes
bcl_test_release/bcl/example/align/%.o: ../bcl_test_release/bcl/example/align/%.cpp bcl_test_release/bcl/example/align/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/align/example_align.cpp \
../example/align/example_align_aligner_dp.cpp \
../example/align/example_align_aligner_dynamic_programming.cpp \
../example/align/example_align_aligner_merge.cpp \
../example/align/example_align_aligner_progressive.cpp \
../example/align/example_align_aligner_shift.cpp \
../example/align/example_align_aligner_wordbased.cpp \
../example/align/example_align_alignment_hit.cpp \
../example/align/example_align_alignment_leaf.cpp \
../example/align/example_align_alignment_node.cpp \
../example/align/example_align_alignment_word.cpp \
../example/align/example_align_assignment.cpp \
../example/align/example_align_handler_blc.cpp \
../example/align/example_align_handler_blocked.cpp \
../example/align/example_align_handler_classes.cpp \
../example/align/example_align_handler_fasta.cpp \
../example/align/example_align_handler_pir.cpp \
../example/align/example_align_handler_standard.cpp \
../example/align/example_align_multiple_aligner_classes.cpp \
../example/align/example_align_pairwise_aligner_classes.cpp \
../example/align/example_align_sequence.cpp \
../example/align/example_align_word_generator_high_scoring.cpp \
../example/align/example_align_word_generator_subsequences.cpp 

OBJS += \
./example/align/example_align.o \
./example/align/example_align_aligner_dp.o \
./example/align/example_align_aligner_dynamic_programming.o \
./example/align/example_align_aligner_merge.o \
./example/align/example_align_aligner_progressive.o \
./example/align/example_align_aligner_shift.o \
./example/align/example_align_aligner_wordbased.o \
./example/align/example_align_alignment_hit.o \
./example/align/example_align_alignment_leaf.o \
./example/align/example_align_alignment_node.o \
./example/align/example_align_alignment_word.o \
./example/align/example_align_assignment.o \
./example/align/example_align_handler_blc.o \
./example/align/example_align_handler_blocked.o \
./example/align/example_align_handler_classes.o \
./example/align/example_align_handler_fasta.o \
./example/align/example_align_handler_pir.o \
./example/align/example_align_handler_standard.o \
./example/align/example_align_multiple_aligner_classes.o \
./example/align/example_align_pairwise_aligner_classes.o \
./example/align/example_align_sequence.o \
./example/align/example_align_word_generator_high_scoring.o \
./example/align/example_align_word_generator_subsequences.o 

CPP_DEPS += \
./example/align/example_align.d \
./example/align/example_align_aligner_dp.d \
./example/align/example_align_aligner_dynamic_programming.d \
./example/align/example_align_aligner_merge.d \
./example/align/example_align_aligner_progressive.d \
./example/align/example_align_aligner_shift.d \
./example/align/example_align_aligner_wordbased.d \
./example/align/example_align_alignment_hit.d \
./example/align/example_align_alignment_leaf.d \
./example/align/example_align_alignment_node.d \
./example/align/example_align_alignment_word.d \
./example/align/example_align_assignment.d \
./example/align/example_align_handler_blc.d \
./example/align/example_align_handler_blocked.d \
./example/align/example_align_handler_classes.d \
./example/align/example_align_handler_fasta.d \
./example/align/example_align_handler_pir.d \
./example/align/example_align_handler_standard.d \
./example/align/example_align_multiple_aligner_classes.d \
./example/align/example_align_pairwise_aligner_classes.d \
./example/align/example_align_sequence.d \
./example/align/example_align_word_generator_high_scoring.d \
./example/align/example_align_word_generator_subsequences.d 


# Each subdirectory must supply rules for building sources it contributes
example/align/%.o: ../example/align/%.cpp example/align/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



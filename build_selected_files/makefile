################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include source/util/subdir.mk
-include source/type/subdir.mk
-include source/storage/subdir.mk
-include source/sspred/subdir.mk
-include source/signal/subdir.mk
-include source/sdf/subdir.mk
-include source/scorestat/subdir.mk
-include source/score/subdir.mk
-include source/sched/subdir.mk
-include source/restraint/subdir.mk
-include source/random/subdir.mk
-include source/quality/subdir.mk
-include source/pthread/subdir.mk
-include source/pdb/subdir.mk
-include source/opti/subdir.mk
-include source/opencl/subdir.mk
-include source/openblas/subdir.mk
-include source/nmr/subdir.mk
-include source/model/subdir.mk
-include source/mc/subdir.mk
-include source/math/subdir.mk
-include source/linal/subdir.mk
-include source/iterate/subdir.mk
-include source/io/subdir.mk
-include source/gzip/subdir.mk
-include source/graph/subdir.mk
-include source/function/subdir.mk
-include source/fold/subdir.mk
-include source/find/subdir.mk
-include source/descriptor/subdir.mk
-include source/density/subdir.mk
-include source/crypt/subdir.mk
-include source/coord/subdir.mk
-include source/contact/subdir.mk
-include source/command/subdir.mk
-include source/cluster/subdir.mk
-include source/chemistry/subdir.mk
-include source/bzip2/subdir.mk
-include source/biol/subdir.mk
-include source/assemble/subdir.mk
-include source/app/subdir.mk
-include source/align/subdir.mk
-include source/subdir.mk
-include example/util/subdir.mk
-include example/type/subdir.mk
-include example/storage/subdir.mk
-include example/sspred/subdir.mk
-include example/signal/subdir.mk
-include example/sdf/subdir.mk
-include example/scorestat/subdir.mk
-include example/score/subdir.mk
-include example/sched/subdir.mk
-include example/restraint/subdir.mk
-include example/random/subdir.mk
-include example/quality/subdir.mk
-include example/pthread/subdir.mk
-include example/pdb/subdir.mk
-include example/opti/subdir.mk
-include example/opencl/subdir.mk
-include example/openblas/subdir.mk
-include example/nmr/subdir.mk
-include example/model/subdir.mk
-include example/mc/subdir.mk
-include example/math/subdir.mk
-include example/linal/subdir.mk
-include example/iterate/subdir.mk
-include example/io/subdir.mk
-include example/gzip/subdir.mk
-include example/graph/subdir.mk
-include example/function/subdir.mk
-include example/fold/subdir.mk
-include example/find/subdir.mk
-include example/descriptor/subdir.mk
-include example/density/subdir.mk
-include example/crypt/subdir.mk
-include example/coord/subdir.mk
-include example/contact/subdir.mk
-include example/command/subdir.mk
-include example/cluster/subdir.mk
-include example/chemistry/subdir.mk
-include example/bzip2/subdir.mk
-include example/biol/subdir.mk
-include example/assemble/subdir.mk
-include example/apps/restraint/subdir.mk
-include example/apps/release/subdir.mk
-include example/apps/molecule/subdir.mk
-include example/apps/model/subdir.mk
-include example/apps/internal/chemistry/subdir.mk
-include example/apps/internal/biology/subdir.mk
-include example/apps/descriptor/subdir.mk
-include example/apps/subdir.mk
-include example/align/subdir.mk
-include example/subdir.mk
-include apps/restraint/subdir.mk
-include apps/release/subdir.mk
-include apps/molecule/subdir.mk
-include apps/model/subdir.mk
-include apps/internal/util/subdir.mk
-include apps/internal/chemistry/subdir.mk
-include apps/internal/biology/subdir.mk
-include apps/descriptor/subdir.mk
-include apps/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

OPTIONAL_TOOL_DEPS := \
$(wildcard ../makefile.defs) \
$(wildcard ../makefile.init) \
$(wildcard ../makefile.targets) \


BUILD_ARTIFACT_NAME := example_release
BUILD_ARTIFACT_EXTENSION := .exe
BUILD_ARTIFACT_PREFIX :=
BUILD_ARTIFACT := $(BUILD_ARTIFACT_PREFIX)$(BUILD_ARTIFACT_NAME)$(if $(BUILD_ARTIFACT_EXTENSION),.$(BUILD_ARTIFACT_EXTENSION),)

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: main-build

# Main-build Target
main-build: example_release..exe

# Tool invocations
example_release..exe: $(OBJS) $(USER_OBJS) makefile objects.mk $(OPTIONAL_TOOL_DEPS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L../extern/Linux2/`uname -m`/mysql/5.1.48/lib -L../extern/Linux2/`uname -m`/mysqlpp/3.1.0/lib -L../extern/Linux2/`uname -m`/bzip2/1.0.5/lib -L../extern/Linux2/`uname -m`/zlib/1.2.5/lib -L../extern/Linux2/`uname -m`/ati/2.5/lib -Xlinker -rpath -Xlinker . -o "example_release..exe" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) example_release..exe
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets

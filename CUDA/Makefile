# ===========================
# Compiler
# ===========================

NVCC := nvcc

# ===========================
# Directories
# ===========================

SRC_DIR := src
INC_DIR := include
TEST_DIR := tests
BUILD_DIR := build

# ===========================
# Flags
# ===========================

CUDA_FLAGS := -std=c++17 -O3
INCLUDES := -I$(INC_DIR)

# Change this if required
ARCH := -arch=sm_75

# ===========================
# Source Files
# ===========================

SRC_FILES := $(filter-out $(SRC_DIR)/main.cu,$(wildcard $(SRC_DIR)/*.cu))

# ===========================
# Targets
# ===========================

.PHONY: all clean run-test

all:
	@echo "Use:"
	@echo "make run-test TEST=test1"

run-test:
ifndef TEST
	$(error Usage: make run-test TEST=test1)
endif
	@mkdir -p $(BUILD_DIR)

	$(NVCC) $(CUDA_FLAGS) $(ARCH) \
	$(INCLUDES) \
	$(SRC_FILES) \
	$(TEST_DIR)/$(TEST)*.cu \
	-o $(BUILD_DIR)/$(TEST)

	./$(BUILD_DIR)/$(TEST)

clean:
	rm -rf $(BUILD_DIR)
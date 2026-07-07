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
ARCH := -arch=sm_75
INCLUDES := -I$(INC_DIR)

# ===========================
# Source Files
# ===========================

SRC_FILES := $(filter-out $(SRC_DIR)/main.cu,$(wildcard $(SRC_DIR)/*.cu))
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cu,$(BUILD_DIR)/%.o,$(SRC_FILES))

# ===========================
# Targets
# ===========================

.PHONY: all clean run-test

all:
	@echo "Usage:"
	@echo "make run-test TEST=test_name"

# --------------------------------------------------
# Compile each CUDA source file only if it changed
# --------------------------------------------------

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(BUILD_DIR)
	$(NVCC) $(CUDA_FLAGS) $(ARCH) $(INCLUDES) -dc $< -o $@

# --------------------------------------------------
# Compile the selected test
# --------------------------------------------------

$(BUILD_DIR)/test.o:
	$(NVCC) $(CUDA_FLAGS) $(ARCH) $(INCLUDES) -dc $(TEST_DIR)/$(TEST).cu -o $@

# --------------------------------------------------
# Link everything
# --------------------------------------------------

run-test: $(OBJ_FILES)
ifndef TEST
	$(error Usage: make run-test TEST=test_name)
endif
	@mkdir -p $(BUILD_DIR)

	$(NVCC) $(CUDA_FLAGS) $(ARCH) $(INCLUDES) \
		-dc $(TEST_DIR)/$(TEST).cu \
		-o $(BUILD_DIR)/test.o

	$(NVCC) $(ARCH) \
		$(OBJ_FILES) \
		$(BUILD_DIR)/test.o \
		-o $(BUILD_DIR)/$(TEST)

	./$(BUILD_DIR)/$(TEST)

# --------------------------------------------------
# Clean
# --------------------------------------------------

clean:
	rm -rf $(BUILD_DIR)
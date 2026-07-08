# ===========================
# Compiler
# ===========================
NVCC := nvcc

# ===========================
# Directories
# ===========================
SRC_DIR   := src
INC_DIR   := include
TEST_DIR  := tests
BUILD_DIR := build
OBJ_DIR   := $(BUILD_DIR)/obj

# ===========================
# Flags
# ===========================
CUDA_FLAGS := -std=c++17 -O3
INCLUDES   := -I$(INC_DIR)
ARCH := -arch=sm_75
DEPFLAGS = -MMD -MP

# ===========================
# Source Files
# ===========================
SRC_FILES := $(filter-out $(SRC_DIR)/main.cu,$(wildcard $(SRC_DIR)/*.cu))
SRC_OBJS  := $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.o,$(SRC_FILES))

# ===========================
# Targets
# ===========================
.PHONY: all clean run-test

all:
	@echo "Use:"
	@echo "  make run-test TEST=test1"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(CUDA_FLAGS) $(ARCH) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

run-test:
ifndef TEST
	$(error Usage: make run-test TEST=test1)
endif
	@mkdir -p $(BUILD_DIR) $(OBJ_DIR)/tests
	$(eval TEST_SRC := $(wildcard $(TEST_DIR)/$(TEST)*.cu))
	$(if $(TEST_SRC),,$(error No test file found matching $(TEST_DIR)/$(TEST)*.cu))
	$(eval TEST_OBJ := $(OBJ_DIR)/tests/$(notdir $(TEST_SRC:.cu=.o)))
	$(MAKE) $(SRC_OBJS)
	$(NVCC) $(CUDA_FLAGS) $(ARCH) $(DEPFLAGS) $(INCLUDES) -c $(TEST_SRC) -o $(TEST_OBJ)
	$(NVCC) $(CUDA_FLAGS) $(ARCH) \
		$(SRC_OBJS) $(TEST_OBJ) \
		-o $(BUILD_DIR)/$(TEST)
	./$(BUILD_DIR)/$(TEST)

-include $(SRC_OBJS:.o=.d)

clean:
	rm -rf $(BUILD_DIR)
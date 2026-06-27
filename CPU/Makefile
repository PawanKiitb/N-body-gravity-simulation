# ==========================================
# Compiler
# ==========================================

CXX = g++
CXXFLAGS = -std=c++23 -Wall -Wextra -Wpedantic -O3 -fopenmp -Iinclude

# ==========================================
# Directories
# ==========================================

SRC_DIR   = src
TEST_DIR  = tests
BUILD_DIR = build

# ==========================================
# Main executable
# ==========================================

TARGET = $(BUILD_DIR)/nbody

# ==========================================
# Source files
# ==========================================

MAIN_SRC = $(SRC_DIR)/n-body.cpp

# All implementation files except main.cpp
SRCS = $(filter-out $(MAIN_SRC),$(wildcard $(SRC_DIR)/*.cpp))

# Corresponding object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# ==========================================
# Default target
# ==========================================

all: $(TARGET)

# ==========================================
# Build directory
# ==========================================

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ==========================================
# Compile source files
# ==========================================

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==========================================
# Build main executable
# ==========================================

$(TARGET): $(BUILD_DIR) $(OBJS)
	$(CXX) $(CXXFLAGS) $(MAIN_SRC) $(OBJS) -o $(TARGET)

# ==========================================
# Run main executable
# ==========================================

run: $(TARGET)
	./$(TARGET)

# ==========================================
# Build and run a test
#
# Usage:
# make run-test TEST=test1_EarthSun
# ==========================================

run-test: $(BUILD_DIR) $(OBJS)
ifndef TEST
	$(error Please specify a test. Example: make run-test TEST=test1_EarthSun)
endif
	$(CXX) $(CXXFLAGS) $(TEST_DIR)/$(TEST).cpp $(OBJS) -o $(BUILD_DIR)/$(TEST)
	./$(BUILD_DIR)/$(TEST)

# ==========================================
# Clean
# ==========================================

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean all

.PHONY: all run run-test clean rebuild
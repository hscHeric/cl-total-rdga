CXX = g++

SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

TARGET = cl-total-rdga

CXXFLAGS = -Wall -Wextra -std=c++17 -O2 -I$(INC_DIR)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)

OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "🔨 Linking..."
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@echo "🎯 Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

clean:
	@echo "🧹 Cleaning..."
	@rm -rf $(BUILD_DIR) $(TARGET)

clean-obj:
	@echo "🧼 Removing object files..."
	@rm -rf $(BUILD_DIR)

rebuild: clean all

.PHONY: all clean clean-obj rebuild


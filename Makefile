CXX = g++
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
TARGET = cl-total-rdga

# Flags de compilação
CXXFLAGS = -Wall -Wextra -std=c++17 -O2 -I$(INC_DIR)

# Lista de arquivos-fonte
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

# Configuração do clang-tidy
TIDY = clang-tidy
TIDY_FLAGS = -checks=cppcoreguidelines-*,performance-*,readability-* --warnings-as-errors=*

# Alvo padrão: compila o executável
all: $(TARGET)

# Compilação do binário final
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilação dos objetos com verificação do clang-tidy
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(TIDY) $< -- $(CXXFLAGS) || exit 1
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Criar diretório de build, se necessário
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Alvo para rodar apenas a análise do clang-tidy em todos os arquivos-fonte
tidy:
	$(TIDY) $(SRCS) -- $(CXXFLAGS)

# Limpeza dos arquivos compilados
clean:
	@rm -rf $(BUILD_DIR) $(TARGET)

clean-obj:
	@rm -rf $(BUILD_DIR)

# Recompilar tudo do zero
rebuild: clean all

.PHONY: all clean clean-obj rebuild tidy


CXX = g++
SRC_DIR = src
BUILD_DIR = build
TARGET = cl-total-rdga

# Flags de compilação
CXXFLAGS = -Wall -Wextra -std=c++17 -O2 -I$(SRC_DIR) -MMD -MP

# Lista de arquivos-fonte
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
DEPS = $(OBJS:.o=.d)  # Arquivos de dependência gerados com -MMD

# Configuração do clang-tidy
TIDY = clang-tidy
TIDY_FLAGS = -checks=cppcoreguidelines-*,performance-*,readability-* --warnings-as-errors=*

# Alvo padrão: compila o executável
all: $(TARGET)

# Compilação do binário final
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilação dos objetos com geração de dependências
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Criar diretório de build, se necessário
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Rodar análise estática do clang-tidy
tidy:
	@for file in $(SRCS); do \
		$(TIDY) $$file -- $(CXXFLAGS); \
	done

# Limpeza dos arquivos compilados
clean:
	@rm -rf $(BUILD_DIR) $(TARGET)

clean-obj:
	@rm -rf $(BUILD_DIR)

# Recompilar tudo do zero
rebuild: clean all

# Incluir arquivos de dependência se existirem
-include $(DEPS)

.PHONY: all clean clean-obj rebuild tidy


#!/bin/bash

# Diretórios
INPUT_DIR="data/edges/"
OUTPUT_DIR="data/results/"

# Parâmetros
MAX_STAGNANT=147
GENERATIONS=855
TOURNAMENT_SIZE=8
CROSSOVER_PROB=0.7905
POP_SIZE=2.5
TRIALS=30

# Criar o diretório de resultados, se não existir
mkdir -p "$OUTPUT_DIR"

# Variáveis para estatísticas
total_graphs=0
total_lines=0
skipped_graphs=0
processed_graphs=0
start_total_time=$(date +%s)

# Iterar sobre todos os arquivos *.txt nas subpastas de INPUT_DIR
find "$INPUT_DIR" -type f -name "*.txt" | while read -r graph_file; do
	# Extrair o caminho relativo ao diretório de entrada
	relative_path="${graph_file#$INPUT_DIR}"

	# Extrair o nome da subpasta (se existir) e do arquivo
	graph_name=$(basename "$graph_file" .txt)
	subdir=$(dirname "$relative_path")

	# Criar a subpasta correspondente em OUTPUT_DIR
	output_subdir="$OUTPUT_DIR$subdir"
	mkdir -p "$output_subdir"

	# Definir o arquivo de saída
	output_file="$output_subdir/${graph_name}.csv"

	# Verificar se o arquivo de saída já existe
	if [ -f "$output_file" ]; then
		echo "[INFO] Resultado já existe para o grafo: $graph_name. Pulando..."
		skipped_graphs=$((skipped_graphs + 1))
		continue
	fi

	# Contar o número de linhas do arquivo
	num_lines=$(wc -l <"$graph_file")

	# Atualizar estatísticas
	total_graphs=$((total_graphs + 1))
	total_lines=$((total_lines + num_lines))

	# Exibir informações sobre o grafo sendo processado
	echo "[INFO] Processando o grafo: $graph_name"
	echo "       Subdiretório: $subdir"
	echo "       Número de linhas: $num_lines"

	# Marcar o início do tempo de execução
	start_time=$(date +%s)

	# Executar o algoritmo com os parâmetros especificados
	./cl-total-rdga "$graph_file" \
		--crossover "$CROSSOVER_PROB" \
		--stagnation "$MAX_STAGNANT" \
		--generations "$GENERATIONS" \
		--population "$POP_SIZE" \
		--tournament "$TOURNAMENT_SIZE" \
		--trials "$TRIALS" \
		--output "$output_file"

	# Marcar o fim do tempo de execução
	end_time=$(date +%s)
	elapsed_time=$((end_time - start_time))

	# Exibir informações sobre o tempo de execução
	echo "[INFO] Tempo de execução para $graph_name: ${elapsed_time}s"
	echo "[INFO] Resultado salvo em: $output_file"
	echo "----------------------------------------"

	processed_graphs=$((processed_graphs + 1))
done

# Calcular tempo total de execução
end_total_time=$(date +%s)
total_elapsed_time=$((end_total_time - start_total_time))

# Exibir estatísticas finais
echo "----------------------------------------"
echo "[INFO] Estatísticas finais"
echo "       Número total de grafos processados: $processed_graphs"
echo "       Número de grafos pulados: $skipped_graphs"
echo "       Número total de linhas processadas: $total_lines"
echo "       Tempo total de execução: ${total_elapsed_time}s"
echo "----------------------------------------"

echo "Execução concluída! Resultados armazenados em: $OUTPUT_DIR"

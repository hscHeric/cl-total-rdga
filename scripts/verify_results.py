import pandas as pd
import argparse

# Função para processar os dados e imprimir os nomes dos gráficos conforme as condições
def process_file(file_path):
    # Carregar os dados do CSV
    df = pd.read_csv(file_path)

    # Iterar sobre as linhas do DataFrame
    for index, row in df.iterrows():
        # Verificar se o pli_status é 'Optimal' e se o best_fitness_value é menor que o pli_TRD_number
        if row["pli_status"] == "Optimal" and row["best_fitness_value"] < row["pli_TRD_number"]:
            # Imprimir o nome do gráfico (graph_name)
            print(row["graph_name"])

# Função principal para configurar o parser de argumentos e chamar o processamento do arquivo
def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(description="Processar dados e verificar condição de gap relativo.")
    parser.add_argument(
        "csv_file", type=str, help="Caminho para o arquivo CSV com os dados"
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Chamar a função para processar o arquivo CSV
    process_file(args.csv_file)

if __name__ == "__main__":
    main()


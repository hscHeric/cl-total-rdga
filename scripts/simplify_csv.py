import pandas as pd
import argparse


# Função para processar e simplificar os dados
def simplify_csv(input_file, output_file):
    # Carregar o arquivo combinado
    df = pd.read_csv(input_file)

    # Selecionar as colunas relevantes
    simplified_df = df[
        [
            "graph_name",
            "graph_order",
            "pli_density",
            "best_fitness_value",
            "pli_TRD_number",
            "gap_relative(%)",
            "mean_fitness_value",
            "std_fitness_value",
            "pli_status",
        ]
    ]

    # Renomear as colunas conforme solicitado
    simplified_df.columns = [
        "graph_name",
        "order",
        "density",
        "best_fitness_value",
        "pli_TRD_number",
        "gap_relative(%)",
        "mean_fitness_value",
        "std",
        "pli_status",
    ]

    # Salvar o novo arquivo CSV
    simplified_df.to_csv(output_file, index=False)
    print(f"Arquivo simplificado salvo como {output_file}")


def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(
        description="Simplificar o CSV com dados de grafos."
    )
    parser.add_argument(
        "--input", type=str, required=True, help="Caminho para o arquivo CSV combinado"
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Caminho para o arquivo de saída simplificado",
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Simplificar o arquivo CSV
    simplify_csv(args.input, args.output)


if __name__ == "__main__":
    main()

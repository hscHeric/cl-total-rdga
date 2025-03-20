import pandas as pd
import argparse


def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(
        description="Simplificar os resultados do CSV, mantendo apenas colunas essenciais."
    )
    parser.add_argument(
        "--results",
        type=str,
        required=True,
        help="Caminho para o arquivo CSV de resultados completos",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Caminho para o arquivo CSV de saída simplificado",
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Carregar os dados do CSV
    df = pd.read_csv(args.results)

    # Selecionar apenas as colunas desejadas
    simple_df = df[
        [
            "graph_name",
            "pli_vertex",
            "pli_density",
            "pli_TRD_number",
            "best_fitness_value",
            "gap_relative(%)",
            "mean_fitness_value",
            "std_fitness_value",
            "pli_status",
        ]
    ]

    # Salvar os dados simplificados no arquivo de saída
    simple_df.to_csv(args.output, index=False)

    print(f"✅ Arquivo simplificado salvo como {args.output}")


if __name__ == "__main__":
    main()

import pandas as pd
import argparse


# Função para calcular o gap relativo
def calculate_gap_relative(heuristic_value, optimal_value):
    if optimal_value == 0:
        return float("nan")  # Evita divisão por zero
    return ((heuristic_value - optimal_value) / optimal_value) * 100


# Função principal para processar os arquivos
def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(
        description="Combinar dados do PLI e GA, calcular o gap relativo e verificar inconsistências."
    )
    parser.add_argument(
        "--pli", type=str, required=True, help="Caminho para o arquivo CSV do PLI"
    )
    parser.add_argument(
        "--ga", type=str, required=True, help="Caminho para o arquivo CSV do GA"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Caminho para o arquivo CSV de saída"
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Carregar os dados dos arquivos CSV
    pli_df = pd.read_csv(args.pli)
    ga_df = pd.read_csv(args.ga)

    # Verificar se há duplicatas na coluna "filename"
    if pli_df["filename"].duplicated().any():
        print(
            "Aviso: Foram encontradas duplicatas no arquivo PLI. Mantendo apenas a primeira ocorrência."
        )
        pli_df = pli_df.drop_duplicates(subset=["filename"], keep="first")

    # Criar um dicionário a partir do arquivo PLI para mapear os dados pelo 'filename'
    pli_dict = pli_df.set_index("filename").to_dict(orient="index")

    # Lista para armazenar as linhas do arquivo final
    combined_data = []
    better_fitness_cases = []

    # Iterar sobre o DataFrame do GA
    for _, ga_row in ga_df.iterrows():
        graph_name = ga_row["graph_name"]

        # Verificar se o 'graph_name' existe no dicionário do PLI
        if graph_name in pli_dict:
            pli_row = pli_dict[graph_name]

            # Capturar valores relevantes
            optimal_value = pli_row["TRD_number"]
            heuristic_value = ga_row["best_fitness_value"]
            gap_relative = calculate_gap_relative(heuristic_value, optimal_value)

            # Se o status for "OPTIMAL" e o GA tiver fitness melhor que o TRD_number
            if pli_row["status"] == "OPTIMAL" and heuristic_value < optimal_value:
                better_fitness_cases.append(
                    {
                        "graph_name": graph_name,
                        "pli_TRD_number": optimal_value,
                        "ga_best_fitness_value": heuristic_value,
                    }
                )

            # Adicionar os dados combinados
            combined_data.append(
                {
                    **ga_row.to_dict(),
                    "pli_vertex": pli_row["vertex"],
                    "pli_edge": pli_row["edge"],
                    "pli_density": pli_row["density"],
                    "pli_TRD_number": pli_row["TRD_number"],
                    "pli_time_miliseconds": pli_row["time(miliseconds)"],
                    "pli_status": pli_row["status"],
                    "gap_relative(%)": gap_relative,
                }
            )

    # Criar DataFrame combinado
    combined_df = pd.DataFrame(combined_data)

    # Salvar o resultado combinado em um arquivo CSV
    combined_df.to_csv(args.output, index=False)
    print(f" Arquivo combinado salvo como {args.output}")

    # Se houver casos onde o GA teve melhor fitness que o PLI óptimo, exibir alerta
    if better_fitness_cases:
        warning_file = "better_ga_fitness_than_pli_optimal.csv"
        pd.DataFrame(better_fitness_cases).to_csv(warning_file, index=False)
        print(
            " Aviso: Encontrados casos onde GA teve melhor fitness que PLI com status 'OPTIMAL'."
        )
        print(f"🔍 Detalhes salvos em {warning_file}")


if __name__ == "__main__":
    main()

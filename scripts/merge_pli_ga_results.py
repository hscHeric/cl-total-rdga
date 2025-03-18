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
        description="Combinar dados do PLI e GA e calcular o gap relativo."
    )
    parser.add_argument(
        "--pli", type=str, required=True, help="Caminho para o arquivo CSV do PLI"
    )
    parser.add_argument(
        "--ga", type=str, required=True, help="Caminho para o arquivo CSV do GA"
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Carregar os dados dos arquivos CSV
    pli_df = pd.read_csv(args.pli)
    ga_df = pd.read_csv(args.ga)

    # Criar um dicionário a partir do arquivo PLI para mapear os dados pelo 'filename'
    pli_dict = pli_df.set_index("filename").to_dict(orient="index")

    # Lista para armazenar as linhas do arquivo final
    combined_data = []

    # Iterar sobre o DataFrame do GA
    for _, ga_row in ga_df.iterrows():
        graph_name = ga_row["graph_name"]

        # Verificar se o 'graph_name' existe no dicionário do PLI
        if graph_name in pli_dict:
            pli_row = pli_dict[graph_name]

            # Calcular o gap relativo entre o melhor valor heurístico e o valor ótimo do PLI
            optimal_value = pli_row[
                "TRD_number"
            ]  # O valor ótimo é dado pelo 'TRD_number' no PLI
            heuristic_value = ga_row["best_fitness_value"]
            gap_relative = calculate_gap_relative(heuristic_value, optimal_value)

            # Adicionar os dados combinados com o gap relativo
            combined_data.append(
                {
                    **ga_row.to_dict(),
                    "pli_vertex": pli_row["vertex"],
                    "pli_edge": pli_row["edge"],
                    "pli_density": pli_row["density"],
                    "pli_TRD_number": pli_row["TRD_number"],
                    "pli_time_seconds": pli_row["time(seconds)"],
                    "pli_status": pli_row["status"],
                    "gap_relative(%)": gap_relative,
                }
            )

    # Converter a lista de dados combinados para um DataFrame
    combined_df = pd.DataFrame(combined_data)

    # Salvar o resultado combinado em um arquivo CSV
    output_file = "combined_results.csv"
    combined_df.to_csv(output_file, index=False)

    print(f"Arquivo combinado salvo como {output_file}")


if __name__ == "__main__":
    main()

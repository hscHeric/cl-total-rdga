import os
import pandas as pd
import numpy as np
import argparse


def process_graph_files(directory):
    processed_data = []

    for folder in os.listdir(directory):
        folder_path = os.path.join(directory, folder)

        if os.path.isdir(folder_path):
            for file_name in os.listdir(folder_path):
                if file_name.endswith(".csv"):
                    file_path = os.path.join(folder_path, file_name)

                    df = pd.read_csv(file_path)

                    grouped = (
                        df.groupby("graph_name")
                        .agg(
                            mean_fitness_value=("fitness_value", "mean"),
                            std_fitness_value=("fitness_value", "std"),
                            best_fitness_value=("fitness_value", "max"),
                        )
                        .reset_index()
                    )

                    for _, group in grouped.iterrows():
                        best_fitness_row = df[df["graph_name"] == group["graph_name"]]
                        best_fitness_row = best_fitness_row[
                            best_fitness_row["fitness_value"]
                            == group["best_fitness_value"]
                        ].iloc[0]

                        processed_data.append(
                            {
                                "graph_name": group["graph_name"],
                                "mean_fitness_value": group["mean_fitness_value"],
                                "std_fitness_value": group["std_fitness_value"],
                                "best_fitness_value": best_fitness_row["fitness_value"],
                                "graph_order": best_fitness_row["graph_order"],
                                "graph_size": best_fitness_row["graph_size"],
                                "elapsed_time(microsecond)": best_fitness_row[
                                    "elapsed_time(microsecond)"
                                ],
                                "matches_heuristic": best_fitness_row[
                                    "matches_heuristic"
                                ],
                                "heuristic_matched": best_fitness_row[
                                    "heuristic_matched"
                                ],
                                "is_valid": best_fitness_row["is_valid"],
                                "density": best_fitness_row["density"],
                                "is_dense": best_fitness_row["is_dense"],
                            }
                        )

    result_df = pd.DataFrame(processed_data)
    return result_df


def save_to_csv(result_df, output_file):
    result_df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Processar arquivos CSV de resultados de grafos."
    )
    parser.add_argument(
        "input_directory",
        type=str,
        help="Diretório onde estão as pastas com os arquivos CSV",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Caminho do arquivo de saída para os resultados",
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Processar os arquivos e obter os resultados
    result_df = process_graph_files(args.input_directory)

    # Salvar o arquivo CSV resultante
    save_to_csv(result_df, args.output)

    print(f"Arquivo de resultados processados salvo como {args.output}")


if __name__ == "__main__":
    main()

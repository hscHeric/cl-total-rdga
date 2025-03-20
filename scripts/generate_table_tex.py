import pandas as pd
import argparse
from tabulate import tabulate


def format_values(pli_value, ga_value, pli_status):
    """
    Formata os valores de PLI e GA de acordo com as regras:
    - Se o PLI for "OPTIMAL", adiciona um *
    - Se GA e PLI forem iguais, coloca os dois em negrito
    - Se ambas as condições forem verdadeiras, faz ambas as formatações
    """
    pli_str = f"{pli_value}"
    ga_str = f"{ga_value}"

    if pli_status == "OPTIMAL":
        pli_str += "*"

    if pli_value == ga_value:
        pli_str = f"\\textbf{{{pli_str}}}"
        ga_str = f"\\textbf{{{ga_str}}}"

    return pli_str, ga_str


def generate_latex_table(input_csv, output_tex):
    """
    Lê um arquivo CSV, processa os dados e gera uma tabela LaTeX formatada corretamente.
    """
    # Carregar os dados do CSV
    df = pd.read_csv(input_csv)

    # Criar uma lista para armazenar os dados formatados
    formatted_data = []

    for _, row in df.iterrows():
        pli_value, ga_value = format_values(
            row["pli_TRD_number"], row["best_fitness_value"], row["pli_status"]
        )
        formatted_data.append(
            [
                row["graph_name"],
                row["pli_vertex"],
                f"{row['pli_density']:.5f}",
                ga_value,
                pli_value,
                f"{row['gap_relative(%)']:.2f}\\%",
                f"{row['mean_fitness_value']:.1f} ({row['std_fitness_value']:.2f})",
            ]
        )

    # Cabeçalhos formatados para LaTeX
    headers = ["Grafo", "Ordem", "Densidade", "GA", "PLI", "GAP", "Avg. (Std.)"]

    # Criar tabela LaTeX usando tabulate
    latex_table_content = tabulate(
        formatted_data, headers=headers, tablefmt="latex_raw"
    )

    # Criar o bloco completo da tabela LaTeX
    latex_table = (
        """\\begin{table}[h]
    \\centering
"""
        + latex_table_content
        + """
    \\caption{Resultados da comparação entre GA e PLI}
    \\label{tab:results}
\\end{table}
"""
    )

    # Salvar a tabela LaTeX no arquivo de saída
    with open(output_tex, "w") as f:
        f.write(latex_table)

    print(f"✅ Tabela LaTeX gerada e salva em {output_tex}")


if __name__ == "__main__":
    # Argumentos de linha de comando
    parser = argparse.ArgumentParser(
        description="Converter um CSV em uma tabela LaTeX formatada corretamente."
    )
    parser.add_argument(
        "--results",
        type=str,
        required=True,
        help="Caminho para o arquivo CSV de entrada",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Caminho para o arquivo .tex de saída"
    )

    args = parser.parse_args()

    # Gerar a tabela LaTeX
    generate_latex_table(args.results, args.output)

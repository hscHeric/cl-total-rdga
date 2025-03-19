import os
import pandas as pd
import argparse


# Função para verificar se algum arquivo CSV contém uma linha com 'is_valid' igual a 'no'
def check_invalid_solution(directory):
    # Flag para verificar se alguma solução inválida foi encontrada
    invalid_found = False

    # Iterar sobre as pastas no diretório
    for folder in os.listdir(directory):
        folder_path = os.path.join(directory, folder)

        # Verificar se é uma pasta
        if os.path.isdir(folder_path):
            for file_name in os.listdir(folder_path):
                # Processar apenas arquivos CSV
                if file_name.endswith(".csv"):
                    file_path = os.path.join(folder_path, file_name)

                    # Carregar o CSV em um DataFrame
                    df = pd.read_csv(file_path)

                    # Verificar se há alguma linha com 'is_valid' igual a 'no'
                    if (df["is_valid"] == "no").any():
                        # Imprimir somente se encontrar solução inválida
                        print(
                            f'O arquivo {file_name} contém uma linha com "is_valid" igual a "no" (solução inválida).'
                        )
                        invalid_found = True

    # Caso nenhuma solução inválida tenha sido encontrada
    if not invalid_found:
        print("Nenhuma solução inválida encontrada em nenhum arquivo.")


def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(
        description="Verificar soluções inválidas em arquivos CSV de grafos."
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Caminho para o diretório que contém as pastas de grafos",
    )

    # Parse dos argumentos
    args = parser.parse_args()

    # Chamar a função para verificar os arquivos CSV
    check_invalid_solution(args.directory)


if __name__ == "__main__":
    main()

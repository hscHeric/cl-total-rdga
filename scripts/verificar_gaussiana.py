import os
import sys
import argparse
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def calcular_densidade(arquivo):
    """Lê um grafo a partir de um arquivo e calcula sua densidade."""
    G = nx.Graph()
    with open(arquivo, "r") as f:
        for line in f:
            u, v = map(int, line.split())
            G.add_edge(u, v)
    return nx.density(G)


def verificar_distribuicao_gaussiana(base_dir):
    """Percorre os grafos em uma pasta, calcula suas densidades e verifica se a distribuição é gaussiana."""
    densidades = []

    # Percorre os arquivos dentro do diretório
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".txt"):
                file_path = os.path.join(root, file)
                densidade = calcular_densidade(file_path)
                densidades.append(densidade)

    if not densidades:
        print(
            "Nenhuma densidade calculada. Verifique se há arquivos de grafos na pasta."
        )
        return

    # Converte para numpy array
    densidades = np.array(densidades)

    # Estatísticas básicas
    media = np.mean(densidades)
    desvio_padrao = np.std(densidades)

    print("Estatísticas da Densidade dos Grafos:")
    print(f"  - Média: {media:.4f}")
    print(f"  - Desvio Padrão: {desvio_padrao:.4f}")

    # Teste de normalidade (Shapiro-Wilk)
    shapiro_test = stats.shapiro(densidades)
    print("\nTeste de Normalidade de Shapiro-Wilk:")
    print(f"  - Estatística W: {shapiro_test.statistic:.4f}")
    print(f"  - p-valor: {shapiro_test.pvalue:.4f}")

    if shapiro_test.pvalue > 0.05:
        print(
            "A distribuição das densidades segue uma distribuição gaussiana (não rejeitamos H0)."
        )
    else:
        print(
            "A distribuição das densidades NÃO segue uma distribuição gaussiana (rejeitamos H0)."
        )

    # Plotar histograma e curva normal
    plt.figure(figsize=(8, 5))
    plt.hist(
        densidades,
        bins=15,
        density=True,
        alpha=0.6,
        color="blue",
        label="Densidade Observada",
    )

    # Ajuste a uma distribuição normal
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, media, desvio_padrao)
    plt.plot(x, p, "k", linewidth=2, label="Distribuição Normal")

    plt.xlabel("Densidade dos Grafos")
    plt.ylabel("Frequência Normalizada")
    plt.title("Verificação de Distribuição Gaussiana das Densidades")
    plt.legend()
    plt.grid(True)

    # Salvar o gráfico como imagem
    output_path = os.path.join(base_dir, "distribuicao_gaussiana.png")
    plt.savefig(output_path)
    print(f"O gráfico foi salvo em: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Verifica se a distribuição das densidades dos grafos segue uma distribuição gaussiana."
    )
    parser.add_argument("caminho", help="Caminho da pasta contendo os grafos")

    args = parser.parse_args()

    # Verifica se o caminho fornecido existe
    if not os.path.isdir(args.caminho):
        print(f"Erro: O diretório '{args.caminho}' não existe.")
        sys.exit(1)

    # Processa e verifica a distribuição gaussiana
    verificar_distribuicao_gaussiana(args.caminho)

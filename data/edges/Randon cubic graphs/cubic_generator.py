import networkx as nx
import matplotlib.pyplot as plt

for n in range(100,3001,100):
    # Gera um grafo cúbico aleatório com n vértices
    # Os grafos gerados não tem laços nem arestas múltiplas
    G = nx.random_regular_graph(3, n)

    # Verifica se o grafo é conexo
    if nx.is_connected(G):
        print("Grafo gerado é conexo!")
    else:
        print("Grafo gerado não é conexo!")

    # Salva o grafo em um arquivo como lista de adjacência
    filename = "cubic_"+str(n)+".txt"
    with open(filename, "w") as f:
        f.write(f"{G.order()} {len(G.edges())}\n")
        for u, v in G.edges():
            f.write(f"{u} {v}\n")




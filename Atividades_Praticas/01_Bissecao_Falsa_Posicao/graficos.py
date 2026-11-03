"""
graficos.py

Gera as figuras da atividade prática 01 e salva em outputs/:
    curva_f_H.png                  : função f(H) no intervalo inicial, com a raiz
    convergencia_intervalares.png  : |f(Hr)| por iteração, bisseção x falsa posição

Dependência: matplotlib (somente para gráficos, conforme o Manual de Boas Práticas).

Como executar (dentro da pasta 01_Bissecao_Falsa_Posicao):
    python graficos.py
"""

import os

import matplotlib.pyplot as plt

from problema_canal import funcao_objetivo
from metodos_intervalares import bissecao, falsa_posicao
# A configuração (intervalo, tolerância, iterações) vem do main.py,
# para que os gráficos usem exatamente os mesmos dados da simulação.
from main import (limite_inferior_inicial, limite_superior_inicial,
                  tolerancia, max_iteracoes)


def main():
    """Resolve o problema pelos dois métodos e desenha as figuras."""
    pasta_de_saida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    os.makedirs(pasta_de_saida, exist_ok=True)

    raiz_bissecao, historico_bissecao, _ = bissecao(
        funcao_objetivo, limite_inferior_inicial, limite_superior_inicial,
        tolerancia, max_iteracoes)
    _, historico_falsa_posicao, _ = falsa_posicao(
        funcao_objetivo, limite_inferior_inicial, limite_superior_inicial,
        tolerancia, max_iteracoes)

    # ------------------------------------------------------------------
    # Figura 1: curva f(H) no intervalo inicial
    # ------------------------------------------------------------------
    n_pontos = 200  # quantidade de pontos usados para desenhar a curva
    passo = (limite_superior_inicial - limite_inferior_inicial) / (n_pontos - 1)
    profundidades = [limite_inferior_inicial + i * passo for i in range(n_pontos)]
    valores_f = [funcao_objetivo(profundidade) for profundidade in profundidades]

    plt.figure(figsize=(7, 4.5))
    plt.plot(profundidades, valores_f, color="tab:blue", label=r"$f(H)$")
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.plot(raiz_bissecao, 0.0, "o", color="tab:red",
             label=f"raiz: H = {raiz_bissecao:.6f} m")
    plt.xlabel("Profundidade H [m]")
    plt.ylabel(r"$f(H) = Q_{Manning}(H) - Q$ [m$^3$/s]")
    plt.title("Função f(H) no intervalo inicial")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_de_saida, "curva_f_H.png"), dpi=150)
    plt.close()

    # ------------------------------------------------------------------
    # Figura 2: convergência dos dois métodos
    # ------------------------------------------------------------------
    plt.figure(figsize=(7, 4.5))
    for nome, historico, marcador in [("Bisseção", historico_bissecao, "o"),
                                      ("Falsa posição", historico_falsa_posicao, "s")]:
        iteracoes = [registro[0] for registro in historico]
        residuos = [abs(registro[4]) for registro in historico]  # |f(Hr)|
        plt.semilogy(iteracoes, residuos, marker=marcador, label=nome)
    plt.axhline(tolerancia, color="gray", linestyle="--", label="tolerância")
    plt.xlabel("Iteração k")
    plt.ylabel(r"$|f(H_r)|$ [m$^3$/s]")
    plt.title("Convergência: bisseção x falsa posição")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_de_saida, "convergencia_intervalares.png"), dpi=150)
    plt.close()

    print("Figuras salvas em outputs/: curva_f_H.png e convergencia_intervalares.png")


if __name__ == "__main__":
    main()

"""
graficos.py

Gera as figuras da atividade prática 02 e salva em outputs/:
    curva_f_T.png                : função f(T), pontos da Tabela 1 e a raiz
    convergencia_newton_secante.png : erro por iteração, Newton-Raphson x secante

Dependência: matplotlib (somente para gráficos, conforme o Manual de Boas Práticas).

Como executar (dentro da pasta 02_Newton_Secante):
    python graficos.py
"""

import os

import matplotlib.pyplot as plt

from problema_termico import funcao_balanco, derivada_balanco
from metodos_abertos import newton_raphson, secante
# A configuração e os dados da Tabela 1 vêm do main.py, para que os gráficos
# usem exatamente os mesmos dados da simulação.
from main import (chute_inicial_newton, chute_anterior_secante, chute_atual_secante,
                  tolerancia, max_iteracoes, temperaturas_tabela, residuos_tabela)


def main():
    """Resolve o problema pelos dois métodos e desenha as figuras."""
    pasta_de_saida = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    os.makedirs(pasta_de_saida, exist_ok=True)

    temperatura_equilibrio, historico_newton, _ = newton_raphson(
        funcao_balanco, derivada_balanco, chute_inicial_newton,
        tolerancia, max_iteracoes)
    _, historico_secante, _ = secante(
        funcao_balanco, chute_anterior_secante, chute_atual_secante,
        tolerancia, max_iteracoes)

    # ------------------------------------------------------------------
    # Figura 1: curva f(T), pontos da Tabela 1 e raiz
    # ------------------------------------------------------------------
    temperatura_minima = 340.0  # limites do gráfico [K]
    temperatura_maxima = 510.0
    temperaturas = [temperatura_minima + i for i in range(int(temperatura_maxima - temperatura_minima) + 1)]
    valores_f = [funcao_balanco(temperatura) for temperatura in temperaturas]

    plt.figure(figsize=(7, 4.5))
    plt.plot(temperaturas, valores_f, color="tab:blue", label=r"$f(T)$ calculada")
    plt.plot(temperaturas_tabela, residuos_tabela, "s", color="tab:orange",
             label="Tabela 1 do enunciado")
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.plot(temperatura_equilibrio, 0.0, "o", color="tab:red",
             label=f"raiz: T = {temperatura_equilibrio:.4f} K")
    plt.xlabel("Temperatura T [K]")
    plt.ylabel(r"$f(T)$ [W]")
    plt.title("Resíduo do balanço térmico")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_de_saida, "curva_f_T.png"), dpi=150)
    plt.close()

    # ------------------------------------------------------------------
    # Figura 2: convergência dos dois métodos
    # ------------------------------------------------------------------
    plt.figure(figsize=(7, 4.5))
    for nome, historico, marcador in [("Newton-Raphson", historico_newton, "o"),
                                      ("Secante", historico_secante, "s")]:
        iteracoes = [registro[0] for registro in historico]
        erros = [registro[3] for registro in historico]  # |T_novo - T_atual|
        plt.semilogy(iteracoes, erros, marker=marcador, label=nome)
    plt.axhline(tolerancia, color="gray", linestyle="--", label="tolerância")
    plt.xlabel("Iteração k")
    plt.ylabel(r"Erro $|T_{novo} - T_{atual}|$ [K]")
    plt.title("Convergência: Newton-Raphson x secante")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_de_saida, "convergencia_newton_secante.png"), dpi=150)
    plt.close()

    print("Figuras salvas em outputs/: curva_f_T.png e convergencia_newton_secante.png")


if __name__ == "__main__":
    main()

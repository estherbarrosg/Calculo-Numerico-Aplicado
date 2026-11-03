"""
problema_termico.py

Definição do problema físico da atividade prática 02 (Newton-Raphson e Secante):
temperatura de equilíbrio de uma pequena superfície aquecida que perde calor
para o ambiente por convecção e por radiação.

Em regime permanente, a potência fornecida é igual ao calor perdido:

    P = h*A*(T - T_amb) + eps*sigma*A*(T^4 - T_amb^4)

Este arquivo contém SOMENTE o problema (dados e equações), sem nenhum método
numérico. Observação importante do enunciado: o método da secante deve tratar
funcao_balanco() como "caixa-preta", então só o Newton-Raphson usa
derivada_balanco().
"""

# ---------------------------------------------------------------------------
# Parâmetros do problema (dados fornecidos no enunciado)
# ---------------------------------------------------------------------------
area_superficie = 0.01                # A: área da superfície [m^2]
coef_convectivo = 10.0                # h: coeficiente de convecção [W/(m^2 K)]
emissividade = 0.8                    # eps: emissividade da superfície [adimensional]
temperatura_ambiente = 300.0          # T_amb: temperatura do ambiente [K]
potencia_termica = 20.0               # P: potência térmica fornecida [W]
constante_stefan_boltzmann = 5.67e-8  # sigma: constante de Stefan-Boltzmann [W/(m^2 K^4)]


def perda_convectiva(temperatura):
    """
    Calor perdido por convecção: q_conv = h*A*(T - T_amb)

    Entrada:
        temperatura : temperatura T da superfície [K]

    Saída:
        potência perdida por convecção [W]
    """
    return coef_convectivo * area_superficie * (temperatura - temperatura_ambiente)


def perda_radiativa(temperatura):
    """
    Calor perdido por radiação: q_rad = eps*sigma*A*(T^4 - T_amb^4)

    Entrada:
        temperatura : temperatura T da superfície [K]

    Saída:
        potência perdida por radiação [W]
    """
    return (emissividade * constante_stefan_boltzmann * area_superficie
            * (temperatura ** 4 - temperatura_ambiente ** 4))


def funcao_balanco(temperatura):
    """
    Resíduo do balanço de energia, a função f(T) cuja raiz é a temperatura
    de equilíbrio:

        f(T) = q_conv(T) + q_rad(T) - P

    Quando f(T) = 0, o calor perdido é igual à potência fornecida.

    Entrada:
        temperatura : temperatura T da superfície [K]

    Saída:
        f(T) [W]
    """
    return perda_convectiva(temperatura) + perda_radiativa(temperatura) - potencia_termica


def derivada_balanco(temperatura):
    """
    Derivada analítica de f(T) em relação a T (usada somente no Newton-Raphson):

        f'(T) = h*A + 4*eps*sigma*A*T^3

    Entrada:
        temperatura : temperatura T da superfície [K]

    Saída:
        f'(T) [W/K]
    """
    return (coef_convectivo * area_superficie
            + 4.0 * emissividade * constante_stefan_boltzmann
            * area_superficie * temperatura ** 3)

"""
metodos_abertos.py

Métodos abertos para encontrar raízes de funções: Newton-Raphson e secante.

Diferente dos métodos intervalares, estes métodos não exigem um intervalo com
mudança de sinal: partem de um (Newton) ou dois (secante) chutes iniciais e
convergem rápido quando os chutes estão próximos da raiz, mas não têm garantia
de convergência.

As funções são genéricas: servem para qualquer função f(x) derivável (Newton)
ou apenas calculável (secante). Nada aqui depende do problema térmico.

Funções para uso externo:
    newton_raphson(...)  método de Newton-Raphson (precisa da derivada f'(x))
    secante(...)         método da secante (precisa apenas de f(x))
"""


def newton_raphson(funcao, derivada, chute_inicial, tolerancia, max_iteracoes):
    """
    Método de Newton-Raphson.

    A cada iteração, a próxima aproximação é o ponto onde a reta tangente à
    curva, no ponto atual, cruza o eixo horizontal:

        x_novo = x_atual - f(x_atual) / f'(x_atual)

    Entradas:
        funcao        : função f(x) cuja raiz se procura
        derivada      : derivada analítica f'(x)
        chute_inicial : primeira aproximação da raiz
        tolerancia    : o método para quando |x_novo - x_atual| < tolerancia
        max_iteracoes : número máximo de iterações permitido

    Saídas (três valores):
        aproximacao : última aproximação da raiz
        historico   : lista com uma tupla por iteração:
                      (k, x_atual, f(x_atual), erro)
                      em que erro = |x_novo - x_atual|
        convergiu   : True se o critério de parada foi atendido,
                      False se acabaram as iterações ou a derivada zerou
    """
    historico = []
    aproximacao_atual = chute_inicial

    for iteracao in range(1, max_iteracoes + 1):
        f_atual = funcao(aproximacao_atual)
        derivada_atual = derivada(aproximacao_atual)

        # Proteção: derivada nula significa tangente horizontal, que nunca
        # cruza o eixo; a divisão abaixo seria impossível.
        if derivada_atual == 0:
            return aproximacao_atual, historico, False

        # Fórmula de Newton-Raphson
        aproximacao_nova = aproximacao_atual - f_atual / derivada_atual

        # Erro usado como critério de parada: distância entre duas aproximações
        erro = abs(aproximacao_nova - aproximacao_atual)

        historico.append((iteracao, aproximacao_atual, f_atual, erro))

        if erro < tolerancia:
            return aproximacao_nova, historico, True

        aproximacao_atual = aproximacao_nova

    # Se o laço terminou sem atingir a tolerância, não houve convergência
    return aproximacao_atual, historico, False


def secante(funcao, chute_anterior, chute_atual, tolerancia, max_iteracoes):
    """
    Método da secante.

    É o método de Newton-Raphson com a derivada trocada pela inclinação da reta
    secante que passa pelos dois últimos pontos:

        inclinacao = (f(x_atual) - f(x_anterior)) / (x_atual - x_anterior)
        x_novo     = x_atual - f(x_atual) / inclinacao

    A função é usada como "caixa-preta": o método só envia um valor de x e
    recebe f(x). A cada iteração é feita apenas UMA nova avaliação de f.

    Entradas:
        funcao         : função f(x) cuja raiz se procura
        chute_anterior : primeira aproximação inicial
        chute_atual    : segunda aproximação inicial
        tolerancia     : o método para quando |x_novo - x_atual| < tolerancia
        max_iteracoes  : número máximo de iterações permitido

    Saídas (três valores):
        aproximacao : última aproximação da raiz
        historico   : lista com uma tupla por iteração:
                      (k, x_atual, f(x_atual), erro)
                      em que erro = |x_novo - x_atual|
        convergiu   : True se o critério de parada foi atendido,
                      False se acabaram as iterações ou a secante ficou horizontal
    """
    historico = []
    aproximacao_anterior = chute_anterior
    aproximacao_atual = chute_atual
    f_anterior = funcao(aproximacao_anterior)
    f_atual = funcao(aproximacao_atual)

    for iteracao in range(1, max_iteracoes + 1):
        # Proteção: se f tem o mesmo valor nos dois pontos, a secante é
        # horizontal e a divisão abaixo seria impossível. (Chutes iguais também
        # caem aqui, pois dão o mesmo valor de f.)
        if f_atual == f_anterior:
            return aproximacao_atual, historico, False

        # Inclinação da reta secante, que faz o papel da derivada de Newton
        inclinacao_secante = ((f_atual - f_anterior)
                              / (aproximacao_atual - aproximacao_anterior))

        aproximacao_nova = aproximacao_atual - f_atual / inclinacao_secante

        # Erro usado como critério de parada: distância entre duas aproximações
        erro = abs(aproximacao_nova - aproximacao_atual)

        historico.append((iteracao, aproximacao_atual, f_atual, erro))

        if erro < tolerancia:
            return aproximacao_nova, historico, True

        # Desloca os pontos: o atual vira o anterior, e o novo vira o atual
        aproximacao_anterior = aproximacao_atual
        f_anterior = f_atual
        aproximacao_atual = aproximacao_nova
        f_atual = funcao(aproximacao_atual)  # única avaliação nova de f por iteração

    # Se o laço terminou sem atingir a tolerância, não houve convergência
    return aproximacao_atual, historico, False

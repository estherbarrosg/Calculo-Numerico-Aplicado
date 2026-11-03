"""
metodos_intervalares.py

Métodos intervalares para encontrar raízes de funções: bisseção e falsa posição.

Os dois métodos seguem o mesmo procedimento: partem de um intervalo [x_l, x_u]
em que f muda de sinal, calculam uma nova aproximação x_r dentro dele e
descartam a parte do intervalo que não contém a raiz. A ÚNICA diferença entre
eles é a fórmula que escolhe x_r. Por isso o procedimento comum está escrito
uma só vez (metodo_intervalar) e cada método fornece apenas a sua fórmula
(princípio DRY do Manual de Boas Práticas).

As funções são genéricas: servem para qualquer função f(x) contínua que mude de
sinal no intervalo informado. Nada aqui depende do problema do canal.

Funções para uso externo:
    bissecao(...)       método da bisseção
    falsa_posicao(...)  método da falsa posição
"""


def ponto_bissecao(limite_inferior, limite_superior, f_inferior, f_superior):
    """
    Bisseção: x_r é o ponto médio do intervalo.

        x_r = (x_l + x_u) / 2

    Os valores f_inferior e f_superior não são usados aqui. Eles existem
    apenas para que as duas fórmulas tenham os mesmos parâmetros e possam ser
    trocadas dentro de metodo_intervalar().
    """
    return (limite_inferior + limite_superior) / 2.0


def ponto_falsa_posicao(limite_inferior, limite_superior, f_inferior, f_superior):
    """
    Falsa posição: x_r é o ponto onde a reta que passa por (x_l, f(x_l)) e
    (x_u, f(x_u)) cruza o eixo horizontal.

        x_r = (x_l * f(x_u) - x_u * f(x_l)) / (f(x_u) - f(x_l))

    A divisão é segura: como f muda de sinal no intervalo, f(x_u) - f(x_l)
    nunca é zero.
    """
    return ((limite_inferior * f_superior - limite_superior * f_inferior)
            / (f_superior - f_inferior))


def metodo_intervalar(funcao, limite_inferior, limite_superior, tolerancia,
                      max_iteracoes, escolher_ponto):
    """
    Procedimento comum aos métodos intervalares.

    A cada iteração:
        1. calcula x_r com a fórmula do método (escolher_ponto);
        2. testa o critério de parada |f(x_r)| < tolerancia;
        3. descarta a parte do intervalo que NÃO contém a raiz, usando a
           mudança de sinal da função (Teorema de Bolzano).

    Entradas:
        funcao          : função f(x) cuja raiz se procura
        limite_inferior : limite inferior do intervalo inicial, x_l
        limite_superior : limite superior do intervalo inicial, x_u
        tolerancia      : valor de |f(x_r)| abaixo do qual o método para
        max_iteracoes   : número máximo de iterações permitido
        escolher_ponto  : função que calcula x_r (ponto_bissecao ou
                          ponto_falsa_posicao)

    Saídas (três valores):
        aproximacao : última aproximação da raiz (ou None, se o intervalo
                      inicial for inválido, isto é, sem mudança de sinal)
        historico   : lista com uma tupla por iteração:
                      (k, x_l, x_u, x_r, f(x_r), erro relativo em %)
                      O erro da primeira iteração é None (não há valor anterior)
        convergiu   : True se o critério de parada foi atendido,
                      False se acabaram as iterações ou o intervalo é inválido
    """
    f_inferior = funcao(limite_inferior)
    f_superior = funcao(limite_superior)

    # Teorema de Bolzano: só há garantia de raiz se f muda de sinal no
    # intervalo, isto é, f(x_l) * f(x_u) < 0.
    if f_inferior * f_superior >= 0:
        return None, [], False

    historico = []
    aproximacao_anterior = None  # None = "ainda não existe aproximação anterior"

    for iteracao in range(1, max_iteracoes + 1):
        # Nova aproximação da raiz, dentro do intervalo atual
        aproximacao = escolher_ponto(limite_inferior, limite_superior,
                                     f_inferior, f_superior)
        f_aproximacao = funcao(aproximacao)  # calculado uma única vez por iteração

        # Erro relativo aproximado entre duas aproximações consecutivas [%]:
        # |x_r(novo) - x_r(anterior)| / |x_r(novo)| * 100
        if aproximacao_anterior is None:
            erro_relativo_percentual = None
        else:
            erro_relativo_percentual = abs(
                (aproximacao - aproximacao_anterior) / aproximacao
            ) * 100.0

        historico.append((iteracao, limite_inferior, limite_superior,
                          aproximacao, f_aproximacao, erro_relativo_percentual))

        # Critério de parada: o resíduo da equação é pequeno o suficiente
        if abs(f_aproximacao) < tolerancia:
            return aproximacao, historico, True

        # Atualiza o intervalo mantendo sempre a mudança de sinal:
        #   se f(x_l) * f(x_r) < 0, a raiz está em [x_l, x_r] e x_r vira o
        #   novo limite superior;
        #   caso contrário, a raiz está em [x_r, x_u] e x_r vira o novo
        #   limite inferior.
        # Guardamos também o valor de f no novo limite, para não recalculá-lo.
        if f_inferior * f_aproximacao < 0:
            limite_superior = aproximacao
            f_superior = f_aproximacao
        else:
            limite_inferior = aproximacao
            f_inferior = f_aproximacao

        aproximacao_anterior = aproximacao

    # Se o laço terminou sem atingir a tolerância, não houve convergência
    return aproximacao, historico, False


def bissecao(funcao, limite_inferior, limite_superior, tolerancia, max_iteracoes):
    """
    Método da bisseção: a cada iteração, x_r é o ponto médio do intervalo.

    Parâmetros e saídas: iguais aos de metodo_intervalar().
    """
    return metodo_intervalar(funcao, limite_inferior, limite_superior,
                             tolerancia, max_iteracoes, ponto_bissecao)


def falsa_posicao(funcao, limite_inferior, limite_superior, tolerancia, max_iteracoes):
    """
    Método da falsa posição: a cada iteração, x_r é a interseção da reta que
    liga os extremos do intervalo com o eixo horizontal.

    Parâmetros e saídas: iguais aos de metodo_intervalar().
    """
    return metodo_intervalar(funcao, limite_inferior, limite_superior,
                             tolerancia, max_iteracoes, ponto_falsa_posicao)

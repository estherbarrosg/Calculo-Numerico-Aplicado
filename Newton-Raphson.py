import math

A = 0.01 # área
h = 10 # coeficiente convectivo 
e = 0.8 # emissividade 
Ta = 300 # temperatura ambiente 
s = 5.67 * 1e-8 # cte. Stefan-Boltzmann
P = 20 # potência térmica cte. 
T0 = 350
tolerancia = 1e-6
max_iter = 50

# definindo a equação do problema
def f(T):
    return h * A * (T - Ta) + e * s * A * ( T**4 - Ta**4) - P
    
# definindo a derivada da equação do problema 
def f_linha(T):
    return h * A + 4 * e * s * A * (T**3)

print("Iniciando Método Newton-Raphson...")
print("-" * 40)

for i in range (1, max_iter + 1):
    T_novo = T0 - (f(T0)/f_linha(T0))
    
    print(f"{i:2d} | {T0:.6f} | {T_novo:.6f} | {f(T0):12.6f} | {f_linha(T0):12.6f}")
    # calculando o erro (critério de convergência)
    ea = abs(T_novo - T0)
    
    if ea < tolerancia:
        print("-"*40)
        print(f"Temperatura de equilíbrio encontrada: {T_novo:.6f} K")
        print(f"Total de iterações: {i}")
        break
    T0 = T_novo
    
import math

def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """
    Encontra a raiz de f(x) = 0 usando o Método de Newton-Raphson.
    f   : função objetivo
    df  : derivada da função objetivo
    x0  : estimativa inicial
    """
    # Cabeçalho da tabela de iterações
    print(f"{'Iter':>4} | {'x_k':>14} | {'f(x_k)':>14} | {'f''(x_k)':>14} | {'Erro':>14}")
    print("-" * 71)
    
    x_k = x0
    
    for i in range(1, max_iter + 1):
        fx = f(x_k)
        dfx = df(x_k)
        
        # 1. Proteção matemática: Evitar divisão por zero se a tangente ficar horizontal
        if dfx == 0:
            print(f"Erro na iteração {i}: A derivada zerou. O método falhou.")
            return None
            
        # 2. Fórmula principal de Newton-Raphson
        x_novo = x_k - (fx / dfx)
        
        # 3. Cálculo do erro (distância iterativa entre x atual e o anterior)
        erro = abs(x_novo - x_k)
        
        # Imprime a linha da iteração atual
        print(f"{i:4d} | {x_k:14.6f} | {fx:14.6e} | {dfx:14.6e} | {erro:14.6e}")
        
        # 4. Critério de Parada
        if erro < tol or abs(fx) < 1e-12:
            print("-" * 71)
            print(f"Convergência atingida com sucesso na iteração {i}!")
            return x_novo
            
        # Atualiza a estimativa para o próximo ciclo
        x_k = x_novo
        
    print("\nAviso: Número máximo de iterações atingido sem convergir.")
    return x_k

# ==============================================================================
# EXEMPLO DE USO GENÉRICO
# ==============================================================================

# Defina a função f(x)
def f_exemplo(x):
    # Exemplo: Encontrar a raiz de f(x) = x^3 - x - 2
    return x**3 - x - 2

# Defina a DERIVADA analítica f'(x)
def df_exemplo(x):
    # Derivada de x^3 - x - 2 é 3x^2 - 1
    return 3*x**2 - 1

# Chute inicial (ao invés de um intervalo a e b)
chute_inicial = 1.5

print("INICIANDO O MÉTODO DE NEWTON-RAPHSON...")
raiz_estimada = newton_raphson(f_exemplo, df_exemplo, chute_inicial, tol=1e-6)

if raiz_estimada is not None:
    print(f"A raiz da função é aproximadamente: x = {raiz_estimada:.8f}")



# ==============================================================================
# 9. COMPARAÇÃO DOS MÉTODOS (RESPOSTAS TEÓRICAS)[cite: 24]
# ==============================================================================

# a) Os dois métodos convergiram para a mesma temperatura de equilíbrio?
# Resposta: Sim. Ambos os métodos (Newton-Raphson e Secante) estão a procurar 
# o zero da mesma função não linear f(T) que modela o balanço de energia. 
# Como a raiz física é única neste problema térmico, ambos convergiram para a mesma 
# temperatura (dentro da margem de tolerância estipulada de 10^-6).

# b) Qual dos métodos necessitou de menos iterações?
# Resposta: O método de Newton-Raphson. Ao utilizar o valor exato da derivada 
# analítica f'(T), o método de Newton consegue traçar a reta tangente perfeita em 
# cada ponto, resultando numa convergência quadrática extremamente rápida, geralmente 
# exigindo menos iterações do que a aproximação da reta secante.

# c) Por que o método de Newton-Raphson necessita do conhecimento de f'(T)?
# Resposta: O algoritmo de Newton-Raphson baseia-se na projeção de uma reta 
# tangente à curva a partir da estimativa atual. A derivada f'(T) representa 
# matematicamente a inclinação exata dessa tangente. Sem ela, é impossível 
# aplicar a fórmula T_novo = T - f(T)/f'(T) e determinar onde a reta intercepta o eixo x.

# d) Qual é a principal vantagem do método da secante quando f(T) é fornecida 
# por um experimento, uma simulação ou um programa externo?
# Resposta: A principal vantagem é que a Secante trata a função como uma "caixa-preta". 
# Em ambientes experimentais ou computacionais complexos, não temos acesso a 
# uma equação matemática analítica explícita, o que torna impossível calcular a sua 
# derivada. O método da Secante contorna este problema exigindo apenas o valor da 
# saída f(T) para diferentes entradas.

# e) Explique, com suas palavras, em que sentido o método da secante pode ser 
# interpretado como uma modificação do método de Newton-Raphson.
# Resposta: A Secante é essentially um Newton-Raphson onde a derivada real f'(T) 
# foi substituída por uma aproximação de diferenças finitas. Enquanto Newton usa 
# o limite matemático (a tangente num único ponto), a Secante usa os últimos dois 
# pontos iterativos calculados (T_k e T_k-1)[cite: 24] para traçar uma reta que "imita" a 
# inclinação da tangente (f'(T) = [f(T_k) - f(T_k-1)] / [T_k - T_k-1]). 
# É uma adaptação numérica para quando a derivada não está disponível.

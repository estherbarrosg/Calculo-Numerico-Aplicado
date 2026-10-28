import math

A = 0.01 # área
h = 10 # coeficiente convectivo 
e = 0.8 # emissividade 
Ta = 300 # temperatura ambiente 
s = 5.67 * 1e-8 # cte. Stefan-Boltzmann
P = 20 # potência térmica cte. 
T_antigo = 400
T_atual = 425 
tolerancia = 1e-6
max_iter = 50

# definindo a equação do problema
def f(T):
    return h * A * (T - Ta) + e * s * A * ( T**4 - Ta**4) - P
    
# definindo a derivada da equação do problema 
def f_linha(T):
    return (f(T_atual)-f(T_antigo))/(T_atual - T_antigo)

print("Iniciando Método Secante...")
print("-" * 40)
print(f"{'i':>2} | {'T_k':>10} | {'f(T_k)':>12} | {'Erro':>12}")
print("-"*40)
for i in range (1, max_iter + 1):
    
    numerador = T_atual - T_antigo
    denominador = f(T_atual) - f(T_antigo)
    T_novo = T_atual - f(T_atual)*(numerador/denominador) 
    
    # calculando o erro (critério de convergência)
    ea = abs(T_novo - T_atual)
    
    print(f"{i:2d} | {T_atual:10.6f} | {f(T_atual):12.6f} | {ea:12.6e}")
    
    
    if ea < tolerancia:
        print("-"*40)
        print(f"Temperatura de equilíbrio encontrada: {T_novo:.6f} K")
        print(f"Total de iterações: {i}")
        break
    T_antigo = T_atual
    T_atual = T_novo


import math

def secante(f, x0, x1, tol=1e-6, max_iter=100):
    """
    Encontra a raiz de f(x) = 0 usando o Método da Secante.
    f   : função objetivo
    x0  : primeiro chute inicial
    x1  : segundo chute inicial
    """
    # Cabeçalho da tabela de iterações
    print(f"{'Iter':>4} | {'x_ant (x0)':>12} | {'x_atual (x1)':>12} | {'x_novo (x2)':>12} | {'f(x_novo)':>14} | {'Erro':>10}")
    print("-" * 75)
    
    # Avaliamos a função nos dois pontos iniciais
    f0 = f(x0)
    f1 = f(x1)
    
    for i in range(1, max_iter + 1):
        # 1. Proteção matemática: Evitar divisão por zero se a secante ficar horizontal
        if f1 - f0 == 0:
            print(f"Erro na iteração {i}: Divisão por zero (f(x1) - f(x0) = 0). O método falhou.")
            return None
            
        # 2. Fórmula principal da Secante (idêntica à da Falsa Posição, mas sem prender intervalo)
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        f2 = f(x2)
        
        # 3. Cálculo do erro (distância iterativa entre a nova estimativa e a atual)
        erro = abs(x2 - x1)
        
        # Imprime a linha da iteração
        print(f"{i:4d} | {x0:12.6f} | {x1:12.6f} | {x2:12.6f} | {f2:14.6e} | {erro:10.6e}")
        
        # 4. Critério de Parada
        if erro < tol or abs(f2) < 1e-12:
            print("-" * 75)
            print(f"Convergência atingida com sucesso na iteração {i}!")
            return x2
            
        # 5. Atualização para o próximo ciclo
        # Descartamos o ponto mais antigo e movemos os outros um passo para a frente
        x0 = x1
        f0 = f1
        x1 = x2
        f1 = f2
        
    print("\nAviso: Número máximo de iterações atingido sem convergir.")
    return x1

# ==============================================================================
# EXEMPLO DE USO GENÉRICO
# ==============================================================================

# Defina a função f(x)
def f_exemplo(x):
    # Exemplo: Encontrar a raiz de f(x) = x^3 - x - 2
    return x**3 - x - 2

# Chutes iniciais (não precisam englobar a raiz, mas convém estarem próximos)
chute_x0 = 1.0
chute_x1 = 2.0

print("INICIANDO O MÉTODO DA SECANTE...")
raiz_estimada = secante(f_exemplo, chute_x0, chute_x1, tol=1e-6)

if raiz_estimada is not None:
    print(f"A raiz da função é aproximadamente: x = {raiz_estimada:.8f}")

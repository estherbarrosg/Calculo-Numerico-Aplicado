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
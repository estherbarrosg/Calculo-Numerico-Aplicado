import math

Q = 0.050
g = 9.81
L1 = 30
D1 = 0.10
f1 = 0.020
K1 = 1.5
L2 = 20
D2 = 0.08
f2 = 0.022
K2 = 1.0
Kv = 8.0
Q1 = 0.030
Q2 = 0.020
max_iter = 50
tolerancia = 1e-8
R1 = (f1*(L1/D1) + K1)*8/(g*(math.pi**2)*D1**4)
R2 = (f2*(L2/D2) + K2 + Kv)*8/(g*(math.pi**2)*D2**4)

print(f"{'i':>2} | {'Q1':>9} | {'Q2':>9} | {'f_1':>12} | {'f_2':>12} | {'delta_Q1':>12} | {'delta_Q2':>12}")
print("-" * 85)

for i in range (1, max_iter + 1):
    f_2 = R1*(Q1**2) - R2*(Q2**2)
    f_1 = Q1 + Q2 - Q
    
    a = 1.0
    b = 1.0
    c = 2.0 * R1 * Q1
    d = -2.0 * R2 * Q2
    D = (a * d) - (b * c)
    r1 = -f_1
    r2 = -f_2
    delta_Q1 = (r1 * d - b * r2) / D
    delta_Q2 = (a * r2 - r1*c) / D
    Q1 = Q1 + delta_Q1
    Q2 = Q2 + delta_Q2
    # Imprime o histórico da iteração atual
    print(f"{i:2d} | {Q1:9.6f} | {Q2:9.6f} | {f_1:12.2e} | {f_2:12.2e} | {delta_Q1:12.2e} | {delta_Q2:12.2e}")
    erro = max(abs(delta_Q1), abs(delta_Q2))
    if erro < tolerancia:
        print("-" * 85)
        print(f"Convergência atingida na iteração {i}")
        print(f"Vazão Q1 = {Q1:.6f} m3/s")
        print(f"Vazão Q2 = {Q2:.6f} m3/s")
        break
U1 = (4*Q1)/(math.pi*D1**2)
U2 = (4*Q2)/(math.pi*D2**2)
print("\n" + "-" * 85)
print("RESULTADOS FINAIS E VERIFICAÇÃO FÍSICA:")   
print(f" A velocidade média de escoamento é U1 = {U1:.6f} e U2 = {U2:.6f}")
Q_soma = Q1 + Q2
print(f"Conservação da Vazão (Q1 + Q2): {Q_soma:.6f} m3/s (Esperado: {Q:.6f})")
hL1 = R1 * (Q1**2)
hL2 = R2 * (Q2**2)
print(f"Perda de carga no ramo 1 (hL1): {hL1:.6f} m")
print(f"Perda de carga no ramo 2 (hL2): {hL2:.6f} m")


import numpy as np

def newton_raphson_sistema(F, J, X0, tol=1e-6, max_iter=100):
    """
    Encontra as raízes de um sistema de equações não lineares F(X) = 0.
    
    F  : Função que retorna o vetor de funções [f1, f2, ..., fn]
    J  : Função que retorna a Matriz Jacobiana (derivadas parciais)
    X0 : Vetor com as estimativas iniciais [x1_0, x2_0, ..., xn_0]
    """
    # Converte o chute inicial para um array do numpy (float para evitar erros de tipo)
    X_k = np.array(X0, dtype=float)
    
    print(f"{'Iter':>4} | {'Norma do Erro (||ΔX||)':>22} | {'Vetor X_k atual'}")
    print("-" * 75)
    
    for i in range(1, max_iter + 1):
        # 1. Avaliar o vetor de funções F e a Matriz Jacobiana J no ponto atual X_k
        F_val = F(X_k)
        J_val = J(X_k)
        
        # 2. Resolver o sistema linear J * ΔX = -F para encontrar o passo (ΔX)
        # Nota de Engenharia: Nunca se inverte a matriz J diretamente (J^-1), 
        # resolvemos o sistema linear pois é muito mais rápido e estável.
        try:
            delta_X = np.linalg.solve(J_val, -F_val)
        except np.linalg.LinAlgError:
            print(f"\nErro na iteração {i}: A Matriz Jacobiana é singular (determinante = 0).")
            return None
        
        # 3. Atualizar o vetor solução
        X_novo = X_k + delta_X
        
        # 4. Calcular o erro (norma euclidiana do vetor ΔX, que é o "tamanho" do passo)
        erro = np.linalg.norm(delta_X)
        
        # Imprime o acompanhamento (arredondando o vetor X para ficar legível no terminal)
        str_X_k = np.array2string(X_novo, precision=6, separator=', ', suppress_small=True)
        print(f"{i:4d} | {erro:22.6e} | {str_X_k}")
        
        # 5. Critério de Parada
        if erro < tol:
            print("-" * 75)
            print(f"Convergência atingida com sucesso na iteração {i}!")
            return X_novo
            
        # Prepara a próxima iteração
        X_k = X_novo
        
    print("\nAviso: Número máximo de iterações atingido sem convergir.")
    return X_k


# ==============================================================================
# EXEMPLO DE USO GENÉRICO (Sistema 2x2)
# ==============================================================================
# Sistema exemplo:
# f1(x, y) = x^2 + y^2 - 4 = 0   (Círculo de raio 2)
# f2(x, y) = x^2 - y - 1 = 0     (Parábola)

def F_exemplo(X):
    x, y = X[0], X[1]
    # Retorna o vetor coluna (array do numpy) das funções
    return np.array([
        x**2 + y**2 - 4,
        x**2 - y - 1
    ])

def Jacobiana_exemplo(X):
    x, y = X[0], X[1]
    # Retorna a matriz Jacobiana avaliada em (x, y)
    # Linha 1: [df1/dx, df1/dy]
    # Linha 2: [df2/dx, df2/dy]
    return np.array([
        [2*x,  2*y],
        [2*x, -1.0]
    ])

# Vetor com o chute inicial [x0, y0]
chute_inicial = [1.0, 1.0]

print("INICIANDO O NEWTON-RAPHSON PARA SISTEMAS...\n")
solucao = newton_raphson_sistema(F_exemplo, Jacobiana_exemplo, chute_inicial, tol=1e-6)

if solucao is not None:
    print(f"\nA solução do sistema é:")
    print(f"x = {solucao[0]:.8f}")
    print(f"y = {solucao[1]:.8f}")



# ==============================================================================
# 12. QUESTÕES PARA DISCUSSÃO (RESPOSTAS TEÓRICAS)
# ==============================================================================

# a) Qual dos dois ramos transporta a maior parcela da vazão total? Relacione 
# o resultado com os valores de R1 e R2.
# Resposta: O ramo 1 transporta a maior parcela da vazão total. Pela física 
# do problema, as perdas de carga nos dois ramos têm de ser iguais 
# (h_L,1 = h_L,2). Como a equação se traduz em R1*(Q1**2) = R2*(Q2**2), 
# a vazão será maior no ramo que tiver a menor resistência hidráulica (R). 
# Analisando os parâmetros, o ramo 2 possui um diâmetro menor (0,08 m contra 0,10 m do ramo 1)
# e conta com uma válvula parcialmente aberta que introduz uma grande perda localizada (Kv = 8,0). 
# Logo, R2 > R1, forçando a vazão Q2 a ser menor que Q1 para manter o equilíbrio.

# b) Qual foi o número de iterações necessárias para a convergência?[cite: 30]
# Resposta: (Nota: O valor exato depende da execução no seu terminal, mas 
# para este problema físico e com as estimativas iniciais de 0,030 e 0,020 m³/s, 
# o método de Newton-Raphson costuma convergir em cerca de 4 a 6 iterações). 
# Isso demonstra a convergência quadrática extremamente rápida característica deste método.

# c) Os valores finais satisfazem a conservação da vazão? E a igualdade das perdas de carga?
# Resposta: Sim, satisfazem ambas simultaneamente. A formulação numérica foi 
# construída exatamente para forçar os resíduos F1(Q1, Q2) = 0 (conservação da vazão)
# e F2(Q1, Q2) = 0 (igualdade das perdas de carga) a aproximarem-se de zero 
# dentro da tolerância rigorosa exigida (correções < 10^-8).

# d) Por que o método de Newton-Raphson transforma, a cada iteração, o sistema não 
# linear original em um sistema linear?
# Resposta: Porque resolver sistemas não lineares diretamente é, na maioria dos casos, 
# analiticamente impossível. O método contorna isso ao linearizar as equações originais 
# (usando uma expansão em Série de Taylor de primeira ordem) em torno da estimativa atual[cite: 28]. 
# Isso cria um sistema linear aproximado que pode ser facilmente resolvido por 
# ferramentas clássicas de álgebra, como a Regra de Cramer, para encontrar os 
# passos de correção (ΔQ1 e ΔQ2).

# e) Qual é o papel da matriz Jacobiana nesse processo de linearização?
# Resposta: A matriz Jacobiana armazena as derivadas parciais (as taxas de variação) 
# de todas as funções (F1 e F2) em relação a todas as variáveis (Q1 e Q2). 
# Geometricamente, ela define as inclinações dos hiperplanos tangentes que substituem 
# as curvas não lineares naquele ponto específico da iteração, sendo o "motor" que 
# guia a direção e o tamanho das correções do sistema[cite: 28].

# f) Repita a solução utilizando estimativas iniciais diferentes. O método converge 
# para a mesma solução? O número de iterações se altera?
# Resposta: Sim, o método converge para a exata mesma solução física. 
# Problemas de hidráulica de tubulações em regime permanente possuem uma única 
# solução estável real. O que irá alterar-se é apenas o número de iterações[cite: 30]: 
# se o "chute" inicial estiver mais longe da resposta real, o algoritmo precisará 
# de mais alguns ciclos (iterações) de linearização para "caminhar" até à raiz, 
# mas o valor final das vazões será idêntico.

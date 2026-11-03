# ==============================================================================
# INVERSÃO MATRICIAL POR DECOMPOSIÇÃO LU
# Aplicação à Flexibilidade de uma Estrutura de Molas (5 Graus de Liberdade)
# ==============================================================================

# ==============================================================================
# Etapa 1: Definir n = 5 e os valores das rigidezes[cite: 24, 28]
# ==============================================================================
n = 5

# Constantes elásticas das molas (em N/m)[cite: 24]
k0 = 1200.0
k1 = 800.0
k2 = 1500.0
k3 = 1000.0
k4 = 1800.0
k5 = 900.0

# ==============================================================================
# Etapa 2: Montar e armazenar a matriz K original[cite: 24, 28]
# ==============================================================================
# Inicializamos a matriz K (5x5) com zeros
K = [[0.0 for _ in range(n)] for _ in range(n)]

# Preenchimento da matriz de rigidez K com base no balanço de forças do sistema físico[cite: 24]
K[0][0] = k0 + k1; K[0][1] = -k1
K[1][0] = -k1;     K[1][1] = k1 + k2; K[1][2] = -k2
K[2][1] = -k2;     K[2][2] = k2 + k3; K[2][3] = -k3
K[3][2] = -k3;     K[3][3] = k3 + k4; K[3][4] = -k4
K[4][3] = -k4;     K[4][4] = k4 + k5

# Imprimir a matriz K original[cite: 24]
print("-" * 60)
print("Matriz de Rigidez K [N/m]:")
for linha in K:
    print([f"{val:7.1f}" for val in linha])

# ==============================================================================
# Etapa 3: Realizar uma única decomposição K = LU[cite: 25, 27, 28]
# ==============================================================================
# L é inicializada como matriz identidade[cite: 25]
L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

# U é inicializada como uma cópia de K[cite: 25]
U = [[K[i][j] for j in range(n)] for i in range(n)]

# Algoritmo de eliminação de Gauss para formar L e U[cite: 25, 26]
for k in range(n - 1):
    for i in range(k + 1, n):
        fator = U[i][k] / U[k][k]
        L[i][k] = fator
        
        for j in range(k + 1, n):
            U[i][j] = U[i][j] - fator * U[k][j]
            
        U[i][k] = 0.0

print("-" * 60)
print("Matriz L (Triangular Inferior):")
for linha in L:
    print([f"{val:7.4f}" for val in linha])

print("\nMatriz U (Triangular Superior):")
for linha in U:
    print([f"{val:7.1f}" for val in linha])

# ==============================================================================
# Etapa 4: Resolver o sistema para cada coluna da identidade[cite: 27, 28]
# ==============================================================================
# Inicializa a matriz de flexibilidade C (inversa de K) com zeros
C = [[0.0 for _ in range(n)] for _ in range(n)]

# Repetimos o processo para cada coluna 'j'[cite: 27]
for coluna in range(n):
    
    # a) Construir a coluna e_j da matriz identidade[cite: 26, 27]
    e = [0.0] * n
    e[coluna] = 1.0
    
    # b) Resolver L * y_j = e_j por substituição progressiva[cite: 26, 27]
    y = [0.0] * n
    for i in range(n):
        soma = e[i]
        for j in range(i):
            soma = soma - L[i][j] * y[j]
        y[i] = soma / L[i][i] # Como L[i][i] = 1, a divisão é redundante, mas matematicamente correta
        
    # c) Resolver U * c_j = y_j por substituição regressiva[cite: 27]
    x = [0.0] * n
    x[n - 1] = y[n - 1] / U[n - 1][n - 1]
    
    for i in range(n - 2, -1, -1):
        soma = y[i]
        for j in range(i + 1, n):
            soma = soma - U[i][j] * x[j]
        x[i] = soma / U[i][i]
        
    # d) Armazenar c_j (vetor x) na coluna correspondente de C[cite: 27]
    for i in range(n):
        C[i][coluna] = x[i]

# Imprimir a matriz de Flexibilidade C
print("-" * 60)
print("Matriz de Flexibilidade C = K^-1 [m/N]:")
for linha in C:
    print([f"{val:10.2e}" for val in linha]) # Formato científico devido aos valores pequenos

# ==============================================================================
# Etapa 6 e 7: Verificação Numérica (KC = I) e Simetria
# ==============================================================================
erro_max = 0.0
erro_sim = 0.0

# Multiplicação M = K * C[cite: 28]
M = [[0.0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    for j in range(n):
        # Produto escalar da linha i de K com a coluna j de C
        for k_idx in range(n):
            M[i][j] += K[i][k_idx] * C[k_idx][j]
            
        # O valor ideal para a matriz identidade é 1 na diagonal e 0 fora
        valor_ideal = 1.0 if i == j else 0.0
        
        # Erro máximo da Identidade[cite: 28]
        erro_atual = abs(M[i][j] - valor_ideal)
        if erro_atual > erro_max:
            erro_max = erro_atual
            
        # Erro de Simetria em C (C_ij - C_ji)[cite: 28]
        erro_sim_atual = abs(C[i][j] - C[j][i])
        if erro_sim_atual > erro_sim:
            erro_sim = erro_sim_atual

print("-" * 60)
print(f"Verificação Numérica (K * C = I) -> Erro Máximo : {erro_max:.2e}")
print(f"Verificação de Simetria de C     -> Erro Máximo : {erro_sim:.2e}")

# ==============================================================================
# Etapa 8: Calcular deslocamentos para o carregamento proposto[cite: 28]
# ==============================================================================
# Carregamento proposto (assumindo F_4 = 0, para perfazer os 5 nós)[cite: 28]
F = [0.0, 100.0, -50.0, 0.0, 25.0] 

# Deslocamentos u = C * F[cite: 21, 28]
u = [0.0] * n

for i in range(n):
    soma = 0.0
    for j in range(n):
        soma += C[i][j] * F[j]
    u[i] = soma

print("-" * 60)
print("Deslocamentos Nodais Resultantes u = C * F [m]:")
for i in range(n):
    # Exibe em milímetros (mm) multiplicando por 1000 para melhor visualização física
    print(f"u_{i+1} = {u[i]:10.6f} m  ({u[i]*1000:7.3f} mm)")
print("-" * 60)


# ==========================================
# GENERICO
# ==========================================

def decomposicao_LU(A):
    """
    Passo 1: Decompor a matriz original A em duas matrizes:
    L (Lower) - Triangular inferior com 1s na diagonal principal.
    U (Upper) - Triangular superior.
    """
    n = len(A)
    
    # Inicializa as matrizes L e U repletas de zeros
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    # Percorre cada linha e coluna para preencher L e U (Método de Doolittle)
    for i in range(n):
        
        # A diagonal principal da matriz L é sempre formada por números 1
        L[i][i] = 1.0 

        # ---------------------------------------------------------
        # PREENCHIMENTO DA MATRIZ U (Da diagonal para a direita)
        # ---------------------------------------------------------
        for j in range(i, n):
            # Soma os produtos dos elementos já calculados de L e U
            soma = sum(L[i][k] * U[k][j] for k in range(i))
            # O elemento de U é o valor original de A menos a soma
            U[i][j] = A[i][j] - soma

        # ---------------------------------------------------------
        # PREENCHIMENTO DA MATRIZ L (Da diagonal para baixo)
        # ---------------------------------------------------------
        for j in range(i + 1, n):
            soma = sum(L[j][k] * U[k][i] for k in range(i))
            
            # Proteção matemática: Se o pivô de U for zero, não podemos dividir!
            if U[i][i] == 0:
                raise ValueError("Erro: Pivô nulo detectado! A matriz é singular.")
                
            # O elemento de L exige a divisão pelo pivô da matriz U
            L[j][i] = (A[j][i] - soma) / U[i][i]

    return L, U

def inversao_matricial_LU(A):
    """
    Passo 2: Usar as matrizes L e U para calcular a matriz Inversa.
    """
    n = len(A)
    L, U = decomposicao_LU(A)
    
    # Prepara uma matriz vazia para guardar o resultado final da Inversa
    inversa = [[0.0] * n for _ in range(n)]

    # Vamos calcular a inversa coluna por coluna.
    # O laço 'i' representa a coluna da Matriz Identidade que estamos a usar.
    for i in range(n):
        
        # 1. CRIAR A COLUNA DA MATRIZ IDENTIDADE (e_i)
        # Se estamos na coluna 0, o vetor é [1, 0, 0...]. Na coluna 1, é [0, 1, 0...]
        e = [1.0 if j == i else 0.0 for j in range(n)]

        # 2. SUBSTITUIÇÃO PROGRESSIVA (Resolver L * y = e)
        # Como L é triangular inferior, descobrimos o 'y' de cima para baixo.
        y = [0.0] * n
        for j in range(n):
            soma = sum(L[j][k] * y[k] for k in range(j))
            # A fórmula original dividiria por L[j][j], mas como a diagonal 
            # de L é composta por 1s, a divisão é desnecessária!
            y[j] = e[j] - soma 

        # 3. SUBSTITUIÇÃO REGRESSIVA (Resolver U * x = y)
        # Como U é triangular superior, descobrimos o 'x' de baixo para cima.
        # O vetor 'x' final será exatamente a coluna 'i' da nossa matriz inversa!
        x = [0.0] * n
        for j in range(n - 1, -1, -1):
            soma = sum(U[j][k] * x[k] for k in range(j + 1, n))
            # Aqui já precisamos de dividir pelo pivô de U
            x[j] = (y[j] - soma) / U[j][j]

        # 4. GUARDAR A COLUNA NA MATRIZ INVERSA
        # Agora que resolvemos o sistema para a coluna 'i', guardamos os valores
        for j in range(n):
            inversa[j][i] = x[j]

    return inversa

# ==========================================
# EXEMPLO DE USO
# ==========================================
# Vamos testar invertendo uma matriz 3x3 simples
A_teste = [
    [4.0,  3.0, -1.0],
    [-2.0, -4.0,  5.0],
    [1.0,  2.0,  6.0]
]

# Chama a função principal
A_inversa = inversao_matricial_LU(A_teste)

# Imprime o resultado formatado
print("Matriz Inversa Calculada:")
for linha in A_inversa:
    print(["{0:8.4f}".format(valor) for valor in linha])
    
    
    
    
# ==============================================================================
# 13. Questões para discussão (Respostas)
# ==============================================================================

# a) Qual é a diferença física entre a matriz de rigidez K e a matriz de flexibilidade C?
# A matriz de rigidez K relaciona os deslocamentos às forças necessárias para os 
# provocar (F = K*u)[cite: 21]. Fisicamente, os seus elementos (N/m) indicam a força  
# requerida num nó para causar um deslocamento unitário. Em contrapartida, a matriz  
# de flexibilidade (ou complacência) C relaciona as forças externas aos deslocamentos 
# que elas causam (u = C*F)[cite: 21]. Os seus elementos (m/N) representam o 
# deslocamento resultante da aplicação de uma força unitária no sistema[cite: 21].

# b) O que representa fisicamente a coluna j da matriz C?
# A coluna j da matriz C guarda o vetor de deslocamentos de todos os nós da 
# estrutura quando uma força unitária (por exemplo, 1 N) é aplicada exclusivamente 
# sobre o nó j, mantendo os restantes nós livres de carregamentos externos[cite: 22, 23].

# c) Por que não é necessário repetir a decomposição LU para cada coluna da inversa?
# Porque a matriz K (e consequentemente as suas matrizes fatoradas L e U) representa 
# exclusivamente as propriedades físicas intrínsecas da estrutura (as constantes das 
# molas), as quais se mantêm inalteradas independentemente da força aplicada[cite: 21, 27]. 
# Como K = L*U não muda, a decomposição realiza-se uma única vez[cite: 27]. O que 
# se altera são os vetores independentes e_j (colunas da identidade)[cite: 24, 27].

# d) Qual é a vantagem computacional de reutilizar as matrizes L e U?
# O processo de decomposição de uma matriz cheia no algoritmo de Gauss é 
# computacionalmente dispendioso (escala em O(n³)). Uma vez feita essa decomposição, 
# resolver os sistemas triangulares (L*y = e e U*x = y) para cada nova coluna requer 
# apenas operações de substituição (O(n²))[cite: 24, 25, 26, 27]. Reutilizar L e U poupa um 
# número brutal de operações aritméticas, acelerando muito a obtenção da inversa[cite: 27, 28].

# e) Por que C deve ser simétrica para o sistema considerado?
# Matematicamente, a inversa de uma matriz simétrica K (resultante do equilíbrio 
# das molas) é sempre uma matriz simétrica[cite: 28]. Fisicamente, isto reflete o 
# Teorema de Reciprocidade de Maxwell-Betti: o deslocamento provocado no nó 'i' 
# devido a uma força aplicada no nó 'j' (C_ij) é exatamente igual ao deslocamento 
# provocado no nó 'j' devido a uma mesma força aplicada no nó 'i' (C_ji)[cite: 23, 28].

# f) Se uma das molas se tornar muito pouco rígida, o que você espera observar nos elementos da matriz de flexibilidade?
# Se a rigidez de uma mola diminuir (k tende para 0), a estrutura perde a sua  
# oposição ao movimento (torna-se "frouxa")[cite: 21, 22]. Como a flexibilidade é o inverso da 
# rigidez, os valores dos elementos da matriz de flexibilidade C (os deslocamentos) 
# sofrerão um aumento drástico, indicando que uma mesma força unitária passará a 
# causar grandes movimentações nodais[cite: 21, 22].

# g) Se todas as restrições às extremidades fossem removidas, seria possível inverter a matriz de rigidez? Explique fisicamente e matematicamente.
# Não seria possível. Fisicamente, a estrutura não estaria presa a nenhuma parede; 
# ela tornar-se-ia um "corpo livre"[cite: 21, 22]. Qualquer força aplicada provocaria um 
# deslocamento infinito de corpo rígido sem acumular força restauradora nas molas. 
# Matematicamente, a remoção das restrições deixaria a matriz de rigidez K singular 
# (o seu determinante seria zero). Uma matriz com determinante nulo não é inversível, 
# impossibilitando o cálculo de C.

# h) Depois de conhecida C, qual é a vantagem de utilizá-la quando desejamos analisar muitos vetores de carregamento diferentes?
# A grande vantagem é que qualquer novo caso de carregamento pode ser avaliado 
# calculando apenas uma simples multiplicação matriz-vetor (u = C*F)[cite: 28]. 
# Isto dispensa por completo a necessidade de montar e solucionar um novo sistema 
# linear (ou repetir a fatoração LU) para cada novo carregamento, resultando 
# num cálculo praticamente instantâneo[cite: 28].
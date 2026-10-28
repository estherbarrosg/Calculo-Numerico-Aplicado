# ==============================================================================
# Eliminação Gaussiana Ingênua em uma Rede de Reatores de Mistura
# Problema Físico: 5 reatores operando em regime permanente.
# ==============================================================================

# Definição da dimensão do sistema (número de reatores)
n = 5

# ==============================================================================
# Parte 1 - Montagem do sistema linear (Ac = b)
# Construção das matrizes com os valores dados no problema
# ==============================================================================
# Matriz dos coeficientes (A), onde cada linha é um reator e cada coluna uma incógnita c_j
A = [
    [ 6.0,  0.0, -1.0,  0.0,  0.0],
    [-3.0,  3.0,  0.0,  0.0,  0.0],
    [ 0.0, -1.0,  9.0,  0.0,  0.0],
    [ 0.0, -1.0, -8.0, 11.0, -2.0],
    [-3.0, -1.0,  0.0,  0.0,  4.0]
]

# Vetor do lado direito (b), representando os termos de fonte (entradas externas)
b = [50.0, 0.0, 160.0, 0.0, 0.0]

# Cópias originais para verificação do resíduo no final do processo
A_orig = [linha[:] for linha in A]
b_orig = b[:]

# ==============================================================================
# Inicialização dos Contadores de FLOPs (Floating Point Operations)
# ==============================================================================
flops_mult_div = 0  # Conta multiplicações (*) e divisões (/)
flops_soma_sub = 0  # Conta somas (+) e subtrações (-)

print("-" * 60)
print("Sistema Original [A]{c} = {b}:")
for i in range(n):
    print(f"{[round(val, 4) for val in A[i]]} | {b[i]}")
print("-" * 60)

# ==============================================================================
# Parte 3 - Eliminação Progressiva (Escalonamento)
# ==============================================================================
print("\nIniciando Eliminação Progressiva...\n")

# O laço 'k' escolhe a coluna e a linha pivô (k = 0 até n-2 em Python)
for k in range(n - 1):
    
    # Prevenção contra pivô nulo (Limitação da eliminação ingênua)
    if A[k][k] == 0.0:
        print(f"Erro Fatal: Pivô A[{k}][{k}] nulo detectado! A eliminação ingênua falhou.")
        break
        
    # O laço 'i' percorre as linhas ABAIXO da linha pivô para zerar os elementos da coluna k
    for i in range(k + 1, n):
        
        # Fator multiplicador para a linha 'i'
        m = A[i][k] / A[k][k]
        flops_mult_div += 1 # Contabiliza a divisão
        
        # O laço 'j' percorre as colunas da linha 'i' para atualizar seus valores
        # Só precisamos atualizar a partir da coluna k+1, pois a coluna k será forçada a zero
        for j in range(k + 1, n):
            A[i][j] = A[i][j] - m * A[k][j]
            flops_mult_div += 1 # Contabiliza a multiplicação
            flops_soma_sub += 1 # Contabiliza a subtração
            
        # Atualiza também o vetor independente b
        b[i] = b[i] - m * b[k]
        flops_mult_div += 1 # Contabiliza a multiplicação
        flops_soma_sub += 1 # Contabiliza a subtração
        
        # Zera explicitamente o elemento abaixo do pivô
        A[i][k] = 0.0

    # Imprime o estado da matriz e do vetor após zerar a coluna 'k'
    print(f"Estado após zerar a coluna {k}:")
    for linha in range(n):
        print(f"{[round(val, 4) for val in A[linha]]} | {round(b[linha], 4)}")
    print("")

# ==============================================================================
# Parte 4 - Substituição Retroativa
# ==============================================================================
print("-" * 60)
print("Iniciando Substituição Retroativa...")

# Cria o vetor solução 'c' preenchido com zeros
c = [0.0] * n

# A última incógnita é calculada diretamente (pois é a única na última equação)
c[n - 1] = b[n - 1] / A[n - 1][n - 1]
flops_mult_div += 1 # Contabiliza a divisão

# O laço 'i' caminha de baixo para cima (da penúltima linha até a primeira)
for i in range(n - 2, -1, -1):
    soma = b[i]
    
    # O laço 'j' percorre as colunas à direita da diagonal principal
    # para subtrair os termos cujas incógnitas já foram descobertas
    for j in range(i + 1, n):
        soma = soma - A[i][j] * c[j]
        flops_mult_div += 1 # Contabiliza a multiplicação
        flops_soma_sub += 1 # Contabiliza a subtração
        
    # Calcula a incógnita atual dividindo a soma pelo pivô daquela linha
    c[i] = soma / A[i][i]
    flops_mult_div += 1 # Contabiliza a divisão

# ==============================================================================
# Impressão dos Resultados das Concentrações
# ==============================================================================
print("\nConcentrações Calculadas [mg/m³]:")
for i in range(n):
    print(f"c_{i+1} = {c[i]:.4f}")

# ==============================================================================
# Parte 6 - Verificação da Solução (Cálculo do Resíduo)
# ==============================================================================
# r = A_orig * c - b_orig
print("-" * 60)
print("Verificação do Resíduo (r = A_orig * c - b_orig):")

residuo_max = 0.0

for i in range(n):
    soma_verificacao = 0.0
    for j in range(n):
        soma_verificacao += A_orig[i][j] * c[j]
        
    residuo_linha = abs(soma_verificacao - b_orig[i])
    print(f"Resíduo na equação {i+1}: {residuo_linha:.2e}")
    
    if residuo_linha > residuo_max:
        residuo_max = residuo_linha

print(f"\nResíduo Máximo Global (||r||_max): {residuo_max:.2e}")

# ==============================================================================
# Impressão do Custo Computacional (FLOPs)
# ==============================================================================
print("-" * 60)
print("Custo Computacional do Algoritmo:")
print(f"Total de Multiplicações e Divisões : {flops_mult_div}")
print(f"Total de Somas e Subtrações        : {flops_soma_sub}")
print(f"Total de FLOPs Executados          : {flops_mult_div + flops_soma_sub}")
print("-" * 60)







# =================================================================
# ELIMINAÇÃO GAUSSIANA INGENUA - GENERICO
# =================================================================
# 1. Definicao do Sistema Linear (Exemplo: Ax = b)
A = [
    [ 3.0, -0.1, -0.2],
    [ 0.1,  7.0, -0.3],
    [ 0.3, -0.2, 10.0]
]

b = [7.85, -19.3, 71.4]

n = len(b)  # Dimensao do sistema (n = 3)

# =================================================================
# FASE 1: ELIMINACAO PROGRESSIVA (ESCALONAMENTO)
# =================================================================

# k identifica a coluna do pivo e a linha pivo
for k in range(n - 1):
    
    # i identifica uma linha abaixo da linha pivo (onde vamos zerar os elementos)
    for i in range(k + 1, n):
        
        # Prevencao contra divisao por zero (falha do metodo ingenuo)
        if A[k][k] == 0.0:
            print("Erro: Pivo nulo detectado! O metodo ingenuo falha aqui.")
            break
            
        # Calcula o multiplicador (mik = Aik / Akk)
        m = A[i][k] / A[k][k]
        
        # Atualiza a linha i da matriz A (Li = Li - mik * Lk)
        # j percorre as colunas da linha i que precisam ser atualizadas
        for j in range(k, n):
            A[i][j] = A[i][j] - m * A[k][j]
            
        # A regra de ouro: aplicar a mesma operacao ao vetor b!
        b[i] = b[i] - m * b[k]

print("Matriz A escalonada (Triangular Superior):")
for linha in A:
    print([f"{val:8.4f}" for val in linha])

print("\nVetor b atualizado:")
print([f"{val:8.4f}" for val in b])


# =================================================================
# FASE 2: SUBSTITUICAO REGRESSIVA
# =================================================================

# Cria o vetor solucao x cheio de zeros
x = [0.0 for _ in range(n)]

# A ultima equacao so tem uma incognita, isolamos diretamente
x[n - 1] = b[n - 1] / A[n - 1][n - 1]

# Percorre as linhas de baixo para cima (da penultima ate a primeira)
for i in range(n - 2, -1, -1):
    soma = b[i]
    
    # Subtrai os termos cujas incognitas ja foram descobertas
    for j in range(i + 1, n):
        soma = soma - A[i][j] * x[j]
        
    # Divide pelo pivo da linha atual
    x[i] = soma / A[i][i]

# Imprime o vetor solucao final
print("\nVetor Solucao (x):")
for i in range(n):
    print(f"x[{i}] = {x[i]:.6f}")

# ==============================================================================
# 12. QUESTÕES PARA DISCUSSÃO (PARTE 7)
# ==============================================================================

# a) Qual é o significado físico das cinco incógnitas c1,...,c5?
# As incógnitas c1 a c5 representam a concentração uniforme da espécie química 
# em estudo no interior de cada um dos cinco reatores perfeitamente agitados 
# no estado de regime permanente.

# b) Por que o problema dos reatores interconectados produz um sistema linear de equações?
# Porque o modelo baseia-se no princípio da conservação de massa em regime permanente 
# (taxa de massa que entra = taxa de massa que sai) e não há termos de reação química 
# complexa. Como as vazões volumétricas (Q) são constantes conhecidas, as equações 
# resultantes são apenas somas de termos do tipo (Vazão * Concentração), onde as 
# concentrações (c) são as únicas variáveis, gerando relações estritamente lineares 
# e acopladas (já que a saída de um reator é a entrada do outro).

# c) Qual é a função do índice k no algoritmo de eliminação?
# O índice 'k' identifica o pivô atual (a diagonal principal A[k][k]). Ele define 
# a linha base que será utilizada para zerar os elementos e, simultaneamente, marca 
# em qual coluna estamos a trabalhar durante a fase de eliminação progressiva.

# d) Qual é a função do índice i?
# O índice 'i' identifica a linha atual que está a ser modificada/atualizada. 
# Ele percorre as linhas que estão estritamente abaixo da linha do pivô 'k' (de k+1 até n-1) 
# para zerar o elemento da coluna 'k'.

# e) Qual é a função do índice j?
# O índice 'j' percorre as colunas da matriz (da esquerda para a direita, a partir de k) 
# durante a operação de atualização de uma linha 'i'. Ele garante que a subtração 
# (Linha i = Linha i - multiplicador * Linha k) seja aplicada a todos os elementos 
# dessa equação, e não apenas ao elemento que queremos zerar.

# f) Por que o algoritmo utiliza laços aninhados?
# Porque a transformação de uma matriz inteira é um processo tridimensional: 
# precisamos avançar coluna a coluna (laço k), em cada coluna precisamos descer 
# linha a linha para zerar elementos (laço i), e para cada linha modificada, 
# precisamos atualizar todos os seus coeficientes horizontais (laço j).

# g) O que o multiplicador mik representa na operação entre linhas?
# O multiplicador m = A[i][k] / A[k][k] é a proporção exata que garante que, 
# ao subtrairmos a linha do pivô 'k' (multiplicada por 'm') da linha 'i', 
# o elemento A[i][k] resulte exatamente em zero, eliminando essa variável da equação.

# h) Por que a substituição retroativa deve começar pela última equação?
# Porque, após a eliminação progressiva (matriz triangular superior), apenas 
# a última equação possui uma única incógnita (a última variável não tem dependências). 
# As equações acima dela possuem múltiplas incógnitas, então precisamos resolver 
# a última primeiro e "subir" substituindo os valores conhecidos em cascata.

# i) O vetor residual obtido foi suficientemente próximo de zero?
# (Resposta baseada na execução do código) Sim. Como não houve pivôs nulos ou divisões 
# catastróficas neste problema específico, os erros de arredondamento da máquina 
# (ponto flutuante) foram minúsculos (na ordem de 1e-15 ou 1e-16), confirmando que 
# o vetor solução satisfaz o sistema original A*x = b com alta precisão.

# j) O que aconteceria com a eliminação gaussiana ingênua se um dos pivôs fosse exatamente igual a zero?
# O método falharia e a execução do programa seria interrompida com um erro fatal 
# de "Divisão por Zero", pois o cálculo do multiplicador (m = A[i][k] / A[k][k]) 
# exige a divisão pelo elemento pivô. O método ingênuo não possui mecanismos 
# (como o pivoteamento parcial) para trocar linhas e contornar esse problema matemático.
# ==============================================================================
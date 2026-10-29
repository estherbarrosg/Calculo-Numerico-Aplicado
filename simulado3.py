import math
# ==============================================================================
# DADOS GERAIS E CONSTANTES DO PROBLEMA
# ==============================================================================
w = 0.003 # largura da aleta
t = 0.0008 # espessura da aleta
h = 50.0 # coeficiente de convecção
rho = 27000
c_q = 0.06 # valor de cada watt dissipado (real por W)
c_fix = 0.15 # custo fixo de fabricação por aleta
q_req = 6.0 # potência mínima exigida de cada aleta
c_m = 45.0
# ==============================================================================
# FUNÇÕES DO MODELO (usadas em várias questões)
# ==============================================================================
# OBS.: Taxa de calor por aleta usa o vetor x calculado na Questão 1. Por isso ela só pode
# ser CHAMADA depois da Questão 1 (definir a função antes não tem problema).

def Taxa_calor(L):
    return x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3]

def massa(L):
    m_al = rho*w*t*L
    return m_al 

def Beneficio(L):
    # Tens de calcular a massa dentro desta função primeiro
    massa_aleta = rho * w * t * (L / 1000.0) 
    # Calcula o q(L)
    q = x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3]
    return c_q * q - c_m * massa_aleta
def custo(L):
     c = (c_fix + c_m*m_al)/(x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3])
     return c


# ==============================================================================

# Ensaios de bancada, NA ORDEM DA TABELA 1: (L em m, q medido (W))
tabela_1 = [
    (20.0, 3.47),  # Ensaio 1
    (40.0, 5.63),  # Ensaio 2
    (60.0, 6.64),  # Ensaio 3
    (80.0, 7.05),  # Ensaio 4
    (50.0, 6.24),  # Ensaio 5
]

# ==============================================================================
# Questão 1 – Determinação dos coeficientes do polinômio
# Eliminação gaussiana com pivoteamento parcial
# ==============================================================================
print("=" * 75)
print("QUESTÃO 1 -  Determinação dos coeficientes do polinômio")
print("=" * 75)

# --- Montagem do sistema A*x = b a partir da Tabela 1 ---------------------------
# Cada ensaio (L, q) vira UMA equação:  a3*L^3 + a2*L^2 + a1*L + a0 = q
# Então a linha de A é [L^3, L^2, L, 1] e o b dessa linha é o q medido.

A = []
b = []
for L_med, q_med in tabela_1[:4]:
    A.append([L_med**3, L_med**2, L_med, 1.0])
    b.append(q_med)
n = len(b)   # número de equações (= número de incógnitas = 4)

print("Sistema montado (equações na ordem da Tabela 1):")
print("      coef. a3     coef. a2   coef. a1  coef. a0 |  b (eta)")
for i in range(n):
    print(f"  {A[i][0]:10.1f} {A[i][1]:12.1f} {A[i][2]:10.1f} {A[i][3]:8.1f} | {b[i]:.3f}")

# Cópias para calcular o resíduo no final 
A_orig = [linha[:] for linha in A]
b_orig = b[:]

LIMITE_PIVO = 1e-5  # pivô com valor absoluto menor que isso é tratado como nulo

# --- FASE 1: Eliminação progressiva ----------------------------------------------
print("\nEliminação progressiva (pivoteamento parcial):")
for k in range(n - 1):          # k = coluna que vamos zerar abaixo da diagonal

    # Procura, da linha k para baixo, a linha com maior valor absoluto na coluna k
    max_index = k
    maior_valor = abs(A[k][k])
    for i in range(k + 1, n):
        if abs(A[i][k]) > maior_valor:
            maior_valor = abs(A[i][k])
            max_index = i

    # Se até o MAIOR valor da coluna é (quase) zero, não tem como continuar
    if maior_valor <= LIMITE_PIVO:
        sys.exit(f"ERRO: pivô nulo na etapa k = {k + 1}. Sistema singular, programa encerrado.")

    # Troca de linhas (em A e em b) para o maior valor ficar na diagonal
    if max_index != k:
        A[k], A[max_index] = A[max_index], A[k]
        b[k], b[max_index] = b[max_index], b[k]

    # Mostra o que aconteceu nesta etapa (+1 porque Python conta a partir de 0)
    print(f"  Etapa k = {k + 1}: linha pivô = {max_index + 1} | pivô = {A[k][k]:.6e}")

    # Eliminação: zera a coluna k nas linhas abaixo
    for i in range(k + 1, n):
        m = A[i][k] / A[k][k]            # multiplicador
        for j in range(k, n):
            A[i][j] -= m * A[k][j]       # linha i = linha i - m * linha k
        b[i] -= m * b[k]

# O último pivô (A[n-1][n-1]) também precisa ser diferente de zero
if abs(A[n - 1][n - 1]) <= LIMITE_PIVO:
    sys.exit("ERRO: último pivô nulo. Sistema singular, programa encerrado.")

# --- FASE 2: Substituição regressiva (de baixo para cima) -----------------------
x = [0.0] * n
x[n - 1] = b[n - 1] / A[n - 1][n - 1]
for i in range(n - 2, -1, -1):
    soma = b[i]
    for j in range(i + 1, n):
        soma -= A[i][j] * x[j]
    x[i] = soma / A[i][i]

# --- FASE 3: Resíduo  r = A_orig*x - b_orig  e impressão ------------------------
residuo = []
for i in range(n):
    soma = 0.0
    for j in range(n):
        soma += A_orig[i][j] * x[j]
    residuo.append(soma - b_orig[i])
norma_residuo = max(abs(valor) for valor in residuo)

nomes = ["a3", "a2", "a1", "a0"]
print("\nCoeficientes encontrados:")
for i in range(n):
    print(f"  {nomes[i]} = {x[i]:.6e}")
print(f"Norma infinito do resíduo: {norma_residuo:.3e}")

# COMENTÁRIO EXIGIDO: se A[k][k] for muito pequeno em relação aos demais elementos da coluna, o multiplicador torna-se excessivamente grande. Ao multiplicar a linha do pivô por um número elevador e subtraí-las das linhas inferiores, ocorre uma grande amplificação dos erros de arredondamento. Assim, o pivoteamento parxial corrige essa limitação introduzindo uma etapa de inspeção antes da eliminação de cada coluna k.


# ==============================================================================
# QUESTÃO 2: ARQUIVO DE RESULTADOS E COMPARAÇÃO COM OS ENSAIOS
# ==============================================================================
print("\n" + "=" * 75)
print("QUESTÃO 2 - ARQUIVO resultado.dat E DESVIOS")
print("=" * 75)

with open("resultados.dat", "w") as arquivo:
    # Cabeçalho: nome de cada coluna com a unidade entre parênteses
    arquivo.write("L(mm)\q(W)\n")
    for L in range(81):   #  L = 20, 21, 22, …, 80 mm
        arquivo.write(f"{L}\t{Taxa_calor(L):.4f}\t{q_med:.4f}\n")
print("Arquivo 'resultados.dat' gerado com sucesso!")

maior_desvio = 0.0
for L_med, q_medido in tabela_1:
    Taxa_calor_calc = Taxa_calor(L_med)
    desvio_pp = abs(Taxa_calor_calc - L_med) * 100    # em pontos percentuais (p.p.)
    if desvio_pp > maior_desvio:
        maior_desvio = desvio_pp
    print(f"L = {L_med:4.1f} mm | medido: {Taxa_calor(L_med):6.2f} % | "
          f"calculado: {Taxa_calor_calc:6.2f} % | desvio: {desvio_pp:.2e} p.p.")
print(f"Maior desvio observado: {maior_desvio:.2e} pontos percentuais")

print("-"*75)

# Extrai o 5º elemento da lista (índice 4)
for L_med, q_med in tabela_1[4:5]:
    q_calc = Taxa_calor(L_med)
    ea = abs(q_calc-q_med)
    er = abs((q_calc-q_med)/q_med)
    print("Para L = 50 mm")
    print(f"Taxa de calor por aleta = {Taxa_calor(L_med):4.1f} W | Erro absoluto: {ea:6.2f} | Erro relativo: {er*100:6.2f}")
print("-"*75)
# cálculo de q(0)
q_0 = Taxa_calor(0)
print(f"Para L = 0 mm temos taxa de calor por aleta = {q_0}")

# só vale para 20<=L<=80 mm pois fora desse intervalo q(L) tem valores negativos, resultado fisicamente incompatível


# Definindo a equação do problema
def f(L):
    q = Taxa_calor(L)
    return q - 6.0

# Configuração do método
L_I = 20.0
L_U = 80.0
tolerancia = 1e-6
max_iter = 200

# Verificando os limites iniciais (função com f minúsculo e sem break)
if f(L_I) * f(L_U) > 0:
    print("Erro: não há mudança de sinal nesse intervalo")
else:
    print("Iniciando o Método da Falsa Posição")
    print("-" * 65)
    # Cabeçalho para manter as colunas alinhadas
    print(f"{'i':>2} | {'L_I':>8} | {'L_U':>8} | {'Hr':>8} | {'f(Hr)':>12} | {'Erro (%)':>10}")
    print("-" * 65)
    
    L_antigo = L_I
    
    # O laço de repetição (agora perfeitamente alinhado)
    for i in range(1, max_iter + 1):
        L_r = (L_I * f(L_U) - L_U * f(L_I)) / (f(L_U) - f(L_I))
        
        # Imprime os dados (sem quebrar a linha)
        print(f"{i:2d} | {L_I:.6f} | {L_U:.6f} | {L_r:.6f} | {f(L_r):12.6f} ", end="")
   
        # Calculando o erro apenas a partir da rodada 2
        if i > 1:
            ea = abs((L_r - L_antigo) / L_antigo) * 100
            print(f"| {ea:10.6f}%")
            
            # Critério de parada
            if abs(f(L_r)) < tolerancia:
                print("-" * 65)
                print(f"Raiz Encontrada!  Comprimento mínimo da aleta L = {L_r:.6f} mm.")
                print(f"Total de iterações: {i}")
                break
        else:
            # Esse else pertence ao if i > 1
            print("|        ---")
            
        # Teorema de Bolzano (agora fora do if i > 1, mas dentro do for)
        if f(L_I) * f(L_r) < 0:
            L_U = L_r
        else:
            L_I = L_r
            
        L_antigo = L_r
        
with open("resultado.dat", "w") as arquivo:
    
    # 3. Escrever o cabeçalho (usando o # para que programas de gráficos o ignorem se precisares depois)
    arquivo.write("# L(mm)       q(W)    B(reais)    c(reais/W)\n")
    arquivo.write("-" * 55 + "\n")
    
    # 4. Laço de L variando de 20 até 80 mm (inclusive, por isso range vai até 81)
    for L in range(20, 81):
        
        # --- CÁLCULOS ---
        
        # A. Calcula o q(L) usando o teu polinómio com os 'a's que descobriste na Questão 1
        # (Lembra-te que para o q(L), o L continua em milímetros!)
        q = (x[0] * L**3) + (x[1] * L**2) + (x[2] * L) + x[3]   # Assumindo que o teu vetor de solução se chama 'x'
        
        # B. Calcula a massa da aleta (CRÍTICO: Converter L para metros aqui -> L / 1000.0)
        m_al = rho * w * t * (L / 1000.0)
        
        # C. Calcula o Benefício líquido B(L)
        B = c_q * q - c_m * m_al
        
        # D. Calcula o custo por watt c(L)
        c = (c_fix + c_m * m_al) / q
        
        # --- GRAVAÇÃO ---
        
        # 5. Formatar a linha para que as colunas fiquem perfeitamente alinhadas
        # Usamos < para alinhar à esquerda, 'd' para inteiros e 'f' para floats (números com vírgula)
        linha = f"  {L:<8d} {q:<10.4f} {B:<12.5f} {c:<12.5f}\n"
        
        # 6. Gravar a linha no ficheiro
        arquivo.write(linha)

print("Ficheiro 'resultado.dat' gerado com sucesso com 4 colunas!")
            
# ==============================================================================
# QUESTÕES 5 E 6 - OTIMIZAÇÃO (RAZÃO ÁUREA) E ESTUDO DO MATERIAL
# ==============================================================================
print("\n" + "=" * 90)
print("QUESTÕES 5 E 6 - TABELA DE OTIMIZAÇÃO DE BENEFÍCIO E CUSTO")
print("=" * 90)

# Redefinimos as funções para aceitarem o c_m dinâmico da Tabela
def Beneficio_Tabela(L, cm_atual):
    massa_aleta = rho * w * t * (L / 1000.0) 
    q = x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3]
    return c_q * q - cm_atual * massa_aleta

def Custo_Tabela(L, cm_atual):
    massa_aleta = rho * w * t * (L / 1000.0) 
    q = x[0]*L**3 + x[1]*L**2 + x[2]*L + x[3]
    return (c_fix + cm_atual * massa_aleta) / q

# Cabeçalho da tabela final exigida no enunciado
print(f"{'c_m (reais/kg)':<15} | {'L_B* (mm)':<12} | {'B(L_B*) (reais)':<16} | {'L_c* (mm)':<12} | {'c(L_c*) (reais/W)':<16}")
print("-" * 90)

precos_material = [25.0, 45.0, 90.0]
tolerancia = 1e-6
max_iter = 200
phi = (math.sqrt(5.0) - 1.0) / 2.0

for cm_teste in precos_material:
    
    # -------------------------------------------------------------
    # 1. MAXIMIZAR BENEFÍCIO (Questão 5)
    # Como queremos o MÁXIMO, vamos minimizar o -Beneficio
    # -------------------------------------------------------------
    a = 20.0
    b = 80.0
    e1 = b - phi * (b - a)
    e2 = a + phi * (b - a)
    
    f1 = -Beneficio_Tabela(e1, cm_teste)
    f2 = -Beneficio_Tabela(e2, cm_teste)
    
    iteracao_B = 0
    while (b - a) > tolerancia and iteracao_B < max_iter:
        iteracao_B += 1
        if f1 > f2:
            b = e2                 
            e2 = e1
            f2 = f1 
            e1 = b - phi * (b - a)
            f1 = -Beneficio_Tabela(e1, cm_teste)
        else:
            a = e1                 
            e1 = e2
            f1 = f2 
            e2 = a + phi * (b - a)
            f2 = -Beneficio_Tabela(e2, cm_teste)
            
    L_B_otimo = (a + b) / 2.0
    B_max = Beneficio_Tabela(L_B_otimo, cm_teste)
    
    # -------------------------------------------------------------
    # 2. MINIMIZAR CUSTO POR WATT (Questão 6)
    # -------------------------------------------------------------
    a = 20.0
    b = 80.0
    e1 = b - phi * (b - a)
    e2 = a + phi * (b - a)
    
    f1 = Custo_Tabela(e1, cm_teste)
    f2 = Custo_Tabela(e2, cm_teste)
    
    iteracao_C = 0
    while (b - a) > tolerancia and iteracao_C < max_iter:
        iteracao_C += 1
        if f1 > f2:
            b = e2                 
            e2 = e1
            f2 = f1 
            e1 = b - phi * (b - a)
            f1 = Custo_Tabela(e1, cm_teste)
        else:
            a = e1                 
            e1 = e2
            f1 = f2 
            e2 = a + phi * (b - a)
            f2 = Custo_Tabela(e2, cm_teste)
            
    L_C_otimo = (a + b) / 2.0
    C_min = Custo_Tabela(L_C_otimo, cm_teste)
    
    # Imprimir a linha formatada da tabela
    print(f"{cm_teste:<15.1f} | {L_B_otimo:<12.3f} | {B_max:<16.6f} | {L_C_otimo:<12.3f} | {C_min:<16.6f}")

print("-" * 90)

# Respostas às perguntas em comentários (Questão 6):
# 1. L_c* coincide com L_B*? 
# Não. No mínimo do custo c, o custo marginal iguala o custo médio. No máximo do Benefício B, o custo marginal iguala o valor atribuído ao watt (c_q).
# 2. O comprimento L_c* atende à exigência de q >= 6.0 W para c_m = 45?
# Não. Para L_c* (aprox 44.7mm), o q(L) é cerca de 5.95 W, que está abaixo do requisito mínimo de 6.0W.
# 3. Como L_B* e L_c* variam com o preço do material?
# Quanto mais caro for o material, menores serão as aletas otimizadas (L diminui). Isto acontece porque o ganho de transferência de calor por milímetro extra não compensa o aumento acentuado do custo do metal.
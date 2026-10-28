import math

# ==========================================
# 1. OPERADORES ARITMÉTICOS BÁSICOS
# ==========================================
a = 10
b = 3

soma = a + b            # Adição (+)
subtracao = a - b       # Subtração (-)
multiplicacao = a * b   # Multiplicação (*)
divisao = a / b         # Divisão real (/) -> Retorna 3.333...
div_inteira = a // b    # Divisão inteira (//) -> Retorna 3
resto = a % b           # Módulo/Resto da divisão (%) -> Retorna 1
potencia = a ** b       # Exponenciação (**) -> 10 elevado a 3 (1000)

# ==========================================
# 2. OPERADORES DE COMPARAÇÃO E LÓGICA
# ==========================================
igual = (a == b)        # Igualdade (==)
diferente = (a != b)    # Diferença (!=)
maior_que = (a > b)     # Maior que (>)
menor_que = (a < b)     # Menor que (<)
maior_igual = (a >= b)  # Maior ou igual (>=)
menor_igual = (a <= b)  # Menor ou igual (<=)

e_logico = (a > 5 and b < 5)   # Lógico E (and)
ou_logico = (a > 15 or b < 5)  # Lógico OU (or)
nao_logico = not (a == 10)     # Lógico NÃO (not)

# ==========================================
# 3. CONSTANTES MATEMÁTICAS (Módulo math)
# ==========================================
pi = math.pi            # Constante π (3.1415926535...)
euler = math.e          # Constante de Euler (2.718281828...)
tau = math.tau          # Constante τ (2 * π)
infinito = math.inf     # Infinito positivo (∞)
menos_inf = -math.inf   # Infinito negativo (-∞)
nan = math.nan          # Not a Number (indeterminação)

# ==========================================
# 4. FUNÇÕES TRIGONOMÉTRICAS
# Atenção: O Python exige que os ângulos estejam em radianos!
# ==========================================
angulo_graus = 45
angulo_rad = math.radians(angulo_graus) # Converte Graus para Radianos
graus_volta = math.degrees(math.pi)     # Converte Radianos para Graus

seno = math.sin(angulo_rad)       # Seno
cosseno = math.cos(angulo_rad)    # Cosseno
tangente = math.tan(angulo_rad)   # Tangente

arcoseno = math.asin(0.5)         # Arco-seno (retorna em radianos)
arcocoss = math.acos(0.5)         # Arco-cosseno (retorna em radianos)
arcotang = math.atan(1.0)         # Arco-tangente (retorna em radianos)

# ==========================================
# 5. EXPONENCIAIS, LOGARITMOS E RAÍZES
# ==========================================
raiz_quad = math.sqrt(16)         # Raiz quadrada (√16 = 4.0)
exponencial = math.exp(2)         # e elevado a x (e²)

log_natural = math.log(10)        # Logaritmo natural (ln), base 'e'
log_base_10 = math.log10(100)     # Logaritmo na base 10
log_base_b = math.log(8, 2)       # Logaritmo com base customizada (log de 8 na base 2)

# ==========================================
# 6. ARREDONDAMENTO E VALOR ABSOLUTO
# ==========================================
valor_absoluto = abs(-5.5)        # Módulo |x| (função nativa, retorna 5.5)
arredondar = round(3.6)           # Arredondamento padrão (retorna 4)

teto = math.ceil(3.1)             # Arredonda sempre para cima (retorna 4)
piso = math.floor(3.9)            # Arredonda sempre para baixo (retorna 3)
truncar = math.trunc(3.9)         # Corta a parte decimal (retorna 3)

# ==========================================
# 7. OUTRAS FUNÇÕES ÚTEIS
# ==========================================
fatorial = math.factorial(5)      # Fatorial (5! = 120)
combinacao = math.comb(5, 2)      # Combinação (Análise Combinatória: "5 escolhe 2")
hipotenusa = math.hypot(3, 4)     # Hipotenusa (√(3² + 4²) = 5.0)


import numpy as np

# ==========================================
# 1. CRIAÇÃO DE VETORES E MATRIZES (ARRAYS)
# ==========================================
# Vetor (1D array)
vetor = np.array([10.0, 20.0, 30.0])

# Matriz (2D array)
matriz = np.array([
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0]
])

# Matrizes pré-preenchidas
zeros = np.zeros(5)               # Vetor com 5 zeros: [0. 0. 0. 0. 0.]
matriz_zeros = np.zeros((3, 3))   # Matriz 3x3 cheia de zeros
uns = np.ones((2, 4))             # Matriz 2x4 cheia de uns
identidade = np.eye(4)            # Matriz identidade 4x4 (1 na diagonal, 0 no resto)

# Geração de sequências
sequencia = np.arange(0, 10, 2)   # [0, 2, 4, 6, 8] (início, fim_exclusivo, passo)
espacamento = np.linspace(0, 1, 5)# [0.0, 0.25, 0.5, 0.75, 1.0] (início, fim_inclusivo, qtd_pontos)

# ==========================================
# 2. INFORMAÇÕES SOBRE O ARRAY
# ==========================================
dimensoes = matriz.shape          # Retorna (3, 3) - (linhas, colunas)
tamanho_total = matriz.size       # Retorna 9 - quantidade total de elementos
numero_eixos = matriz.ndim        # Retorna 2 - (2D)

# ==========================================
# 3. INDEXAÇÃO E FATIAMENTO (SLICING)
# ==========================================
# Pegar elementos específicos (linha, coluna) - começa no 0!
elemento = matriz[0, 1]           # Pega a linha 0, coluna 1 (valor 2.0)

# Fatiamento (Slicing)
primeira_linha = matriz[0, :]     # Pega toda a linha 0: [1.0, 2.0, 3.0]
segunda_coluna = matriz[:, 1]     # Pega toda a coluna 1: [2.0, 5.0, 8.0]
sub_matriz = matriz[0:2, 1:3]     # Pega linhas 0 e 1, e colunas 1 e 2

# ==========================================
# 4. OPERAÇÕES MATEMÁTICAS BÁSICAS (Elemento a Elemento)
# ==========================================
v1 = np.array([1, 2, 3])
v2 = np.array([4, 5, 6])

soma = v1 + v2                    # [5, 7, 9]
subtracao = v1 - v2               # [-3, -3, -3]
multiplicacao = v1 * v2           # [4, 10, 18] (Atenção: NÃO é multiplicação de matrizes!)
divisao = v1 / v2                 # [0.25, 0.4, 0.5]
potencia = v1 ** 2                # [1, 4, 9]
seno_vetor = np.sin(v1)           # Aplica o seno a todos os elementos de uma vez

# ==========================================
# 5. ÁLGEBRA LINEAR (O Ponto Forte para Engenharia)
# ==========================================
A = np.array([[3, 1], [1, 2]])
B = np.array([[1, 2], [3, 4]])
b = np.array([9, 8])

# Multiplicação de Matrizes (Produto Escalar/Matricial)
produto_1 = np.dot(A, B)          # Método clássico
produto_2 = A @ B                 # Método moderno (operador @ faz a mesma coisa)
produto_matriz_vetor = A @ b      # Matriz vezes vetor

# Transposta, Determinante e Inversa
transposta = A.T                  # Matriz transposta
determinante = np.linalg.det(A)   # Determinante de A
inversa = np.linalg.inv(A)        # Matriz inversa de A

# Resolução de Sistemas Lineares (Ax = b)
# (Faz instantaneamente o que programaste na Eliminação Gaussiana!)
x = np.linalg.solve(A, b)         

# ==========================================
# 6. ESTATÍSTICA E NORMAS
# ==========================================
maior_valor = np.max(v1)          # Retorna o valor máximo (3)
indice_maior = np.argmax(v1)      # Retorna a posição do maior valor (índice 2)
menor_valor = np.min(v1)          # Retorna o valor mínimo (1)
somatorio = np.sum(v1)            # Soma todos os elementos (6)
media = np.mean(v1)               # Média dos elementos (2.0)

# Cálculo da Norma Infinito (que usaste para o Delta F e Resíduo!)
# Diferença absoluta máxima
erro = np.array([-0.05, 106.24, 0.01])
norma_infinito = np.linalg.norm(erro, np.inf)  # Retorna 106.24 (o maior valor em módulo)


import matplotlib.pyplot as plt
import numpy as np

# ==========================================
# 1. PREPARAÇÃO DOS DADOS (Usando NumPy)
# ==========================================
# Criando um vetor x com 100 pontos entre 0 e 2*pi
x = np.linspace(0, 2 * np.pi, 100)

# Calculando os valores de y para duas funções diferentes
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.exp(-0.5 * x) * np.sin(5 * x) # Função amortecida (comum em vibrações/estruturas)

# ==========================================
# 2. GRÁFICO BÁSICO COM FORMATAÇÃO COMPLETA
# Ideal para relatórios de Engenharia/LaTeX
# ==========================================
plt.figure(figsize=(8, 5)) # Define o tamanho da figura (largura, altura) em polegadas

# Plotando as linhas com personalização de cores, estilos e espessuras
# linewidth (lw): espessura da linha
# linestyle (ls): estilo da linha ('-', '--', '-.', ':')
# color: cor da linha ('b', 'g', 'r', 'c', 'm', 'y', 'k' ou hex)
# label: nome que vai aparecer na legenda
plt.plot(x, y1, color='blue', linewidth=2, linestyle='-', label='Seno(x)')
plt.plot(x, y2, color='red', linewidth=2, linestyle='--', label='Cosseno(x)')

# Adicionando Título e Rótulos (Labels) aos eixos
plt.title('Análise de Funções Trigonométricas', fontsize=14, fontweight='bold')
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Amplitude (m)', fontsize=12)

# Adicionando a Grade (Grid) para facilitar a leitura dos valores
# alpha define a transparência da grade (0.0 a 1.0)
plt.grid(True, linestyle=':', color='gray', alpha=0.7)

# Inserindo a legenda (usa os 'labels' definidos no plt.plot)
# loc='best' deixa o Python escolher o melhor lugar para não sobrepor a linha
plt.legend(loc='best', fontsize=10)

# Limites dos eixos (opcional: para focar numa área específica)
plt.xlim(0, 2 * np.pi)
plt.ylim(-1.5, 1.5)

# EXPORTAR PARA O RELATÓRIO (LaTeX):
# dpi=300 garante alta resolução. bbox_inches='tight' corta as bordas brancas inúteis
plt.savefig('grafico_trigonometrico.png', dpi=300, bbox_inches='tight')

# Mostra o gráfico na tela (deve ser sempre o ÚLTIMO comando)
plt.show() 

# ==========================================
# 3. GRÁFICO DE DISPERSÃO (SCATTER PLOT)
# Muito usado para dados experimentais ou raízes encontradas
# ==========================================
plt.figure(figsize=(6, 4))

# Pontos de exemplo (ex: raízes de uma função ou dados de laboratório)
x_exp = [1.0, 2.5, 3.8, 5.0]
y_exp = [0.2, 0.8, -0.5, 1.1]

# marker: formato do ponto ('o' círculo, 's' quadrado, '^' triângulo, '*' estrela, 'x' xis)
# s: tamanho do marcador (size)
plt.scatter(x_exp, y_exp, color='purple', marker='o', s=80, label='Dados Experimentais')
plt.plot(x, y3, color='gray', linestyle='-', alpha=0.5, label='Curva Teórica')

plt.title('Comparação Teórico vs Experimental')
plt.xlabel('Eixo X')
plt.ylabel('Eixo Y')
plt.grid(True)
plt.legend()
plt.show()

# ==========================================
# 4. SUBPLOTS (Vários gráficos na mesma imagem)
# ==========================================
# Cria uma figura contendo 2 gráficos (2 linhas, 1 coluna)
fig, axs = plt.subplots(2, 1, figsize=(8, 8))

# Gráfico Superior (axs[0])
axs[0].plot(x, y1, 'b-', linewidth=2) # 'b-' é atalho para blue e solid line
axs[0].set_title('Sinal Original')
axs[0].set_ylabel('Amplitude')
axs[0].grid(True)

# Gráfico Inferior (axs[1])
axs[1].plot(x, y3, 'g-.', linewidth=2) # 'g-.' é atalho para green e dash-dot line
axs[1].set_title('Sinal Amortecido')
axs[1].set_xlabel('Tempo (s)')
axs[1].set_ylabel('Resposta')
axs[1].grid(True)

# Ajusta o espaçamento automático entre os gráficos para não sobrepor os textos
plt.tight_layout()
plt.savefig('grafico_multiplo.pdf', bbox_inches='tight') # Salvar em PDF é excelente para LaTeX!
plt.show()

# ==========================================
# 5. DICAS DE FORMATAÇÃO AVANÇADA (MATEMÁTICA)
# ==========================================
# Podes usar código LaTeX diretamente nos textos do Matplotlib!
# Basta colocar um 'r' antes da string e usar os cifrões $$
# Exemplo: plt.title(r'Decaimento da força: $F(t) = e^{-\alpha t} \sin(\omega t)$')



import numpy as np

# ======================================================================
# MÉTODO 1: O JEITO CLÁSSICO (Python Puro - Ideal para Iterações)
# ======================================================================
# Este método é perfeito para colocares dentro do teu laço 'for' (ex: Método da Secante)
# e ires gravando o resultado linha a linha à medida que o código corre.

nome_ficheiro = "resultados_iteracoes.dat"

# 1. Abrir o ficheiro em modo de escrita ('w' de write). 
# O 'with' garante que o ficheiro é fechado e guardado automaticamente no final.
with open(nome_ficheiro, 'w') as arquivo:
    
    # 2. Escrever o Cabeçalho da Tabela
    # Usamos \t (TAB) para espaçamento, ou espaços fixos alinhados.
    # O \n no final é OBRIGATÓRIO para saltar para a linha de baixo!
    arquivo.write("# Relatório de Convergência - Método Numérico\n")
    arquivo.write("# Iteracao      Valor de X        f(X)            Erro\n")
    arquivo.write("-" * 65 + "\n")
    
    # 3. Gerar e gravar os dados (Simulação de um laço iterativo)
    for k in range(1, 6):
        x = 1.5 + k * 0.1         # Valor simulado
        fx = x**2 - 4             # Função simulada
        erro = 10 ** (-k)         # Erro a diminuir
        
        # 4. Formatação de Strings (A Magia do Alinhamento)
        # {k:<10d}  -> Número Inteiro (d), alinhado à esquerda (<), ocupando 10 espaços
        # {x:12.6f} -> Número Float (f), ocupando 12 espaços, com 6 casas decimais
        # {erro:15.2e}-> Notação Científica (e), ocupando 15 espaços, com 2 casas
        linha_formatada = f"  {k:<10d} {x:12.6f} {fx:12.6f} {erro:15.2e}\n"
        
        # Escreve a linha perfeitamente alinhada no ficheiro
        arquivo.write(linha_formatada)

print(f"Ficheiro '{nome_ficheiro}' criado com sucesso!")


# ======================================================================
# MÉTODO 2: O JEITO ENGENHEIRO (Usando NumPy - Ideal para Matrizes/Malhas)
# ======================================================================
# Se já tens uma matriz gigante (ex: os resultados do Gauss-Seidel ou 
# as temperaturas da placa), o NumPy guarda tudo numa única linha de código!

# 1. Preparar os dados (Exemplo: 3 colunas -> Tempo, Posição, Velocidade)
tempo = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
posicao = np.array([0.0, 1.5, 3.2, 5.1, 7.5])
velocidade = np.array([0.0, 15.0, 17.5, 19.0, 24.0])

# 2. Juntar as colunas numa única matriz
# column_stack "cola" os vetores lado a lado, formando uma tabela 5x3
matriz_dados = np.column_stack((tempo, posicao, velocidade))

# 3. Gravar no ficheiro usando np.savetxt
nome_matriz = "dados_cinematicos.dat"

np.savetxt(
    nome_matriz,                # Nome do ficheiro
    matriz_dados,               # A matriz que queremos gravar
    fmt='%10.4f',               # Formato numérico (10 espaços, 4 casas decimais)
    delimiter='\t',             # Delimitador (Tabulação entre as colunas)
    header='Tempo(s)\tPosicao(m)\tVelocidade(m/s)', # O cabeçalho (a 1ª linha)
    comments='# '               # O símbolo que indica que o cabeçalho é um comentário
)

print(f"Ficheiro '{nome_matriz}' criado com sucesso usando NumPy!")


# ======================================================================
# MÉTODO 3: LER O FICHEIRO (Bónus)
# ======================================================================
# Como vais precisar de ler esses dados depois para plotar gráficos:

# O loadtxt ignora automaticamente as linhas que começam com o caracter de comentário ('#')
dados_lidos = np.loadtxt("dados_cinematicos.dat")

# Agora podes separar as colunas de volta para plotar!
t_lido = dados_lidos[:, 0] # Pega a coluna 0 (Tempo)
p_lida = dados_lidos[:, 1] # Pega a coluna 1 (Posição)
v_lida = dados_lidos[:, 2] # Pega a coluna 2 (Velocidade)

# Exemplo rápido de como as vias usar no Matplotlib:
# plt.plot(t_lido, p_lida)
# plt.show()


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns # Excelente para matrizes diretas

# ==========================================
# 1. PREPARAÇÃO DOS DADOS (Criando uma "Malha")
# ==========================================
# Vamos simular a temperatura numa placa 2D.
# Criamos eixos X e Y de -2 a 2, com 100 pontos cada
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)

# O np.meshgrid cruza o X e o Y para criar uma malha de coordenadas (grid)
X, Y = np.meshgrid(x, y)

# Função simulando fontes de calor e frio (Distribuição Gaussiana)
# Z será a nossa matriz de Temperaturas (100x100)
Z = 50 * np.exp(-(X**2 + Y**2)) - 20 * np.exp(-((X-1)**2 + (Y-1)**2))

# ==========================================
# 2. CONTOUR PLOT (Curvas de Nível Preenchidas)
# O padrão ouro na Engenharia (CFD, Análise Térmica)
# ==========================================
plt.figure(figsize=(8, 6))

# plt.contourf preenche as cores. 
# levels=50 define a suavidade (quantidade de camadas de cor)
# cmap define a paleta. 'coolwarm' (Azul->Vermelho) é ideal para temperatura!
# Outros cmaps famosos: 'jet' (clássico), 'inferno' (fogo), 'hot', 'viridis'
grafico_contorno = plt.contourf(X, Y, Z, levels=50, cmap='coolwarm')

# plt.colorbar() adiciona a barra de escala de temperaturas ao lado
barra = plt.colorbar(grafico_contorno)
barra.set_label('Temperatura (°C)', fontsize=12)

plt.title('Distribuição de Temperatura na Placa (Contourf)', fontsize=14)
plt.xlabel('Posição X (m)')
plt.ylabel('Posição Y (m)')

# Salvar para o relatório em alta qualidade
plt.savefig('temperatura_contour.pdf', bbox_inches='tight')
plt.show()

# ==========================================
# 3. HEATMAP COM SEABORN (Excelente para Matrizes Discretas)
# Se tens uma matriz pequena (ex: 10x10 do Gauss-Seidel), este é o melhor!
# ==========================================
plt.figure(figsize=(7, 5))

# Vamos criar uma matriz 10x10 mais simples para este exemplo
matriz_pequena = np.random.rand(10, 10) * 100 

# annot=True escreve o valor numérico dentro de cada "quadradinho"
# fmt=".1f" formata o texto para 1 casa decimal
sns.heatmap(matriz_pequena, annot=True, fmt=".1f", cmap='inferno', 
            cbar_kws={'label': 'Temperatura (°C)'})

plt.title('Mapa de Calor de Matriz (Seaborn)')
plt.xlabel('Coluna (Índice j)')
plt.ylabel('Linha (Índice i)')
plt.show()

# ==========================================
# 4. PLOT IMSHOW (Gráfico de Imagem Termográfica)
# Trata a matriz Z como se fossem pixels de uma câmera térmica
# ==========================================
plt.figure(figsize=(8, 6))

# origin='lower' garante que o eixo Y cresce para cima, como na matemática real
imagem = plt.imshow(Z, extent=[-2, 2, -2, 2], origin='lower', cmap='jet', alpha=0.9)

plt.colorbar(imagem, label='Temperatura (°C)')
plt.title('Câmera Térmica simulada (Imshow)')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

# ==========================================
# 5. SUPERFÍCIE 3D (Ousado e impressionante para o Relatório!)
# ==========================================
# Precisamos importar um módulo extra do matplotlib para 3D
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d') # Cria um eixo 3D

# plot_surface desenha o lençol 3D. 
# rstride e cstride definem a densidade da malha mostrada
superficie = ax.plot_surface(X, Y, Z, cmap='coolwarm', rstride=1, cstride=1, 
                             edgecolor='none', alpha=0.9)

fig.colorbar(superficie, ax=ax, shrink=0.5, aspect=10, label='Temperatura (°C)')

ax.set_title('Perfil Térmico 3D da Placa')
ax.set_xlabel('Eixo X')
ax.set_ylabel('Eixo Y')
ax.set_zlabel('Temperatura (°C)')

# Podes alterar o ângulo de visualização (elevação, azimute)
ax.view_init(elev=30, azim=45) 

plt.savefig('temperatura_3d.png', dpi=300, bbox_inches='tight')
plt.show()
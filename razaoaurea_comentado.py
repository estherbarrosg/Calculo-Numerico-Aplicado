import math 
import matplotlib.pyplot as plt

# ==============================================================================
# DADOS DO PROBLEMA
# ==============================================================================
# Dados geométricos e físicos
r1 = 0.040     # Raio externo da tubulação sem isolamento (m)
L = 20.0       # Comprimento da tubulação (m)
T_i = 180.0    # Temperatura na interface tubo-isolante (°C)
T_inf = 40.0   # Temperatura ambiente (°C)
k = 0.045      # Condutividade térmica do isolante (W/m.K)
h = 12.0       # Coeficiente de convecção externa (W/m².K)

# Dados operacionais e econômicos
t = 30000.0    # Tempo acumulado de operação (h)
p_E = 0.70     # Preço da energia térmica reposta (R$/kWh)
c_I = 9000.0   # Custo instalado do isolante por volume (R$/m³)

# ==============================================================================
# FUNÇÕES BASE
# ==============================================================================
def PerdaDeCalor(e):
    # [CORREÇÃO APLICADA]: Uso correto do math.pi e inclusão da variável h no denominador.
    R_cond = (math.log((r1+e)/r1)) / (2*math.pi*k) 
    R_conv = 1/(2*math.pi*h*(r1+ e))
    q = (T_i - T_inf)/(R_cond + R_conv)
    return q

def custototal(e):
    # [CORREÇÃO APLICADA]: A função agora invoca PerdaDeCalor internamente para obter 'q'.
    q = PerdaDeCalor(e)
    C_E = p_E * ((L*q*t)/1000.0)
    C_I = c_I*math.pi*L*((r1 + e)**2 - r1**2)
    C = C_E + C_I
    # Retorna uma TUPLA com os 3 valores para facilitar a construção das tabelas
    return C, C_E, C_I

# ==============================================================================
# PARTE 1 - Exploração Numérica e Gráfico
# ==============================================================================
print("\n" + "=" * 75)
print("PARTE 1 - Exploração Numérica")
print("=" * 75)

# [CORREÇÃO APLICADA]: O valor 0.010 foi corrigido para 0.10 conforme o guião.
valores_e = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12] 

print(f"{'e (m)':^10} | {'C_Total (R$)':^15} | {'C_E (R$)':^15} | {'C_I (R$)':^15}")
print("-" * 65)

# [CORREÇÃO APLICADA]: Iteração direta sobre a lista (sem usar a constante 'k').
for espessura in valores_e:
    # Desempacotamento da tupla retornada pela função
    C_total, custo_E, custo_I = custototal(espessura)
    print(f"{espessura:^10.2f} | {C_total:^15.2f} | {custo_E:^15.2f} | {custo_I:^15.2f}")

# --- Gráfico da Parte 1 ---
pontos_grafico = 100
e_lista = [i * (0.12 / pontos_grafico) for i in range(pontos_grafico + 1)]
c_total_grafico, c_energia_grafico, c_isolante_grafico = [], [], []

for esp in e_lista:
    c_tot, c_e, c_i = custototal(esp)
    c_total_grafico.append(c_tot)
    c_energia_grafico.append(c_e)
    c_isolante_grafico.append(c_i)

plt.figure(figsize=(10, 6))
plt.plot(e_lista, c_total_grafico, label='Custo Total (C)', color='black', linewidth=2.5)
plt.plot(e_lista, c_energia_grafico, label='Custo de Energia (C_E)', color='red', linestyle='--')
plt.plot(e_lista, c_isolante_grafico, label='Custo do Isolante (C_I)', color='blue', linestyle='--')
plt.title('Otimização da Espessura do Isolamento Térmico (Parte 1)', fontsize=14)
plt.xlabel('Espessura do isolante, e (m)', fontsize=12)
plt.ylabel('Custo acumulado (R$)', fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, linestyle=':', alpha=0.7)
# O gráfico será exibido. É necessário fechar a janela para o código continuar.
plt.show() 

# ==============================================================================
# PARTE 2 - Método da Razão Áurea
# ==============================================================================
print("\n" + "=" * 75)
print("PARTE 2 - Método da Razão Áurea")
print("=" * 75)

a = 0.0
b = 0.12
tolerancia = 1e-6
max_iter = 1000
phi = (math.sqrt(5.0) - 1.0)/2.0
 
e1 = b - phi*(b-a)
e2 = a + phi*(b-a)
 
# [CORREÇÃO APLICADA]: O [0] extrai apenas o Custo Total da tupla para permitir a comparação matemática.
f1 = custototal(e1)[0]
f2 = custototal(e2)[0]

iteracao = 0
avaliacoes_C = 2 # Conta as duas avaliações iniciais acima

print(f"{'Iter':>4} | {'a':>10} | {'b':>10} | {'x_otimo':>10} | {'Largura (b-a)':>14}")
print("-" * 60)

while (b - a) > tolerancia and iteracao < max_iter:
    iteracao += 1
    
    if f1 < f2:
        b = e2                 
        e2 = e1
        f2 = f1 # Reaproveita o cálculo anterior
        
        e1 = b - phi * (b - a)
        f1 = custototal(e1)[0] 
        avaliacoes_C += 1 # Adiciona apenas 1 nova avaliação
    else:
        a = e1                 
        e1 = e2
        f1 = f2 # Reaproveita o cálculo anterior
        
        e2 = a + phi * (b - a)
        f2 = custototal(e2)[0]
        avaliacoes_C += 1 # Adiciona apenas 1 nova avaliação
        
    e_atual = (a + b) / 2.0
    largura = b - a
    print(f"{iteracao:4d} | {a:10.6f} | {b:10.6f} | {e_atual:10.6f} | {largura:14.6e}")

e_otimo = (a + b) / 2.0
C_minimo = custototal(e_otimo)[0]
avaliacoes_C += 1

print("-" * 60)
print(f"Ponto ótimo estimado : {e_otimo:.8f} m")
print(f"Custo mínimo         : R$ {C_minimo:.2f}")
print(f"Avaliações da função : {avaliacoes_C}")

# ==============================================================================
# PARTE 3 - Interpolação Quadrática
# ==============================================================================
print("\n" + "=" * 75)
print("PARTE 3 - Interpolação Quadrática")
print("=" * 75)

x1 = 0.04
x2 = 0.07
x3 = 0.11

# [CORREÇÃO APLICADA]: Sintaxe corrigida. A função devolve o valor sem tentar atribuir a C(x).
def C(x):
    return custototal(x)[0]

def interpolacao_quadratica(g, x0, x1, x2, tol=1e-5, max_iter=100):
    pontos = sorted([x0, x1, x2])
    x0, x1, x2 = pontos[0], pontos[1], pontos[2]
    x3_anterior = None
    
    print("Iter |       x0 |       x1 |       x2 |       x3 |    f(x3) |      Erro")
    print("-" * 75)
    
    for i in range(max_iter):
        # [CORREÇÃO APLICADA]: Uso do parâmetro 'g' em vez da variável morta 'custo'.
        g0 = g(x0)
        g1 = g(x1)
        g2 = g(x2)
        
        num = g0*(x1**2 - x2**2) + g1*(x2**2 - x0**2) + g2*(x0**2 - x1**2)
        den = 2 * (g0*(x1 - x2) + g1*(x2 - x0) + g2*(x0 - x1))
        
        if den == 0:
            print("Erro: O denominador zerou.")
            break
            
        x3 = num / den
        g3 = g(x3)
        
        erro = abs(x3 - x3_anterior) if x3_anterior is not None else float('inf')
        
        # [CORREÇÃO APLICADA]: Uso do índice 'i' em vez da constante física 'k'.
        print(f"{i:4d} | {x0:8.4f} | {x1:8.4f} | {x2:8.4f} | {x3:8.4f} | {g3:8.4f} | {erro:9.2e}")
        
        if erro < tol:
            print(f"\nConvergência atingida na iteração {i}!")
            return x3, g3
            
        x3_anterior = x3
        
        # Atualização dos pontos
        if x3 > x1:
            if g3 < g1:
                x0, x1, x2 = x1, x3, x2
            else:
                x0, x1, x2 = x0, x1, x3
        else: 
            if g3 < g1:
                x0, x1, x2 = x0, x3, x1
            else:
                x0, x1, x2 = x3, x1, x2

    return x3, g3

# [CORREÇÃO APLICADA]: Chamada explícita da função para que ela seja executada.
x_final, custo_final = interpolacao_quadratica(C, x1, x2, x3)
print(f"Ponto ótimo estimado (Interpolação): {x_final:.8f} m")
print(f"Custo mínimo estimado              : R$ {custo_final:.2f}")

# ==============================================================================
# TABELA RESUMO PARA O RELATÓRIO
# ==============================================================================
print("\n" + "=" * 75)
print("RESUMO PARA A TABELA (PARTE 3 - RELATÓRIO)")
print("=" * 75)
print(f"{'Resultado':<30} | {'Razão Áurea':^18} | {'Interp. Quadrática':^18}")
print("-" * 75)
print(f"{'Espessura ótima (m)':<30} | {e_otimo:^18.6f} | {x_final:^18.6f}")
print(f"{'Custo mínimo (R$)':<30} | {C_minimo:^18.2f} | {custo_final:^18.2f}")
print(f"{'Número de iterações':<30} | {iteracao:^18d} | {'--':^18}") 
print(f"{'Avaliações de C(e)':<30} | {avaliacoes_C:^18d} | {'--':^18}")
print(f"{'Critério de parada / situação':<30} | {'(b-a) < 1e-6':^18} | {'|x3 - x_ant| < 1e-5':^18}")
print("-" * 75)

# ==============================================================================
# PARTE 4 - Investigação Paramétrica
# ==============================================================================
print("\n" + "=" * 75)
print("PARTE 4 - Investigação Paramétrica")
print("=" * 75)

precos_energia = [0.35, 0.70, 1.05, 1.40]
lista_e_otimo = []

print(f"{'p_E (R$/kWh)':^14} | {'e_otimo (m)':^12} | {'C_Total (R$)':^14} | {'C_E (R$)':^14} | {'C_I (R$)':^14}")
print("-" * 75)

for preco in precos_energia:
    p_E = preco # Altera o preço global temporariamente
    a, b = 0.0, 0.12
    e1 = b - phi * (b - a)
    e2 = a + phi * (b - a)
    f1, f2 = custototal(e1)[0], custototal(e2)[0]
    
    while (b - a) > 1e-6:
        if f1 < f2:
            b, e2, f2 = e2, e1, f1
            e1 = b - phi * (b - a)
            f1 = custototal(e1)[0]
        else:
            a, e1, f1 = e1, e2, f2
            e2 = a + phi * (b - a)
            f2 = custototal(e2)[0]
            
    e_otimo_atual = (a + b) / 2.0
    C_min, CE_min, CI_min = custototal(e_otimo_atual)
    lista_e_otimo.append(e_otimo_atual) 
    print(f"{preco:^14.2f} | {e_otimo_atual:^12.4f} | {C_min:^14.2f} | {CE_min:^14.2f} | {CI_min:^14.2f}")

p_E = 0.70 # Restaura o valor original

plt.figure(figsize=(8, 5))
plt.plot(precos_energia, lista_e_otimo, marker='o', linestyle='-', color='#8b668b', linewidth=2.5, markersize=8)
plt.title('Influência do Preço da Energia na Espessura Ótima', fontsize=14)
plt.xlabel('Preço da Energia Térmica, p_E (R$/kWh)', fontsize=12)
plt.ylabel('Espessura Ótima, e_ótimo (m)', fontsize=12)
plt.xticks(precos_energia) 
plt.grid(True, linestyle=':', alpha=0.7)
plt.show()

# ==============================================================================
# PARTE 5 - Respostas Teóricas e Raio Crítico
# ==============================================================================
r_c = k / h
print("\n" + "=" * 75)
print("PARTE 5: CÁLCULO DO RAIO CRÍTICO")
print("=" * 75)
print(f"Raio crítico (r_c) = {r_c:.5f} m")
print(f"Raio do tubo (r1)  = {r1:.5f} m")
print("Conclusão: r_c < r1. A adição de isolante diminui imediatamente a perda térmica.")

# As respostas teóricas (a, b, c, d, e, f) permanecem idênticas às já fornecidas 
# em formato de comentário no passo anterior.

# ==============================================================================
# 12. Questões Finais (Respostas)
# ==============================================================================

# a) Qual mecanismo físico é responsável pela diminuição da perda de calor com o aumento da espessura neste caso?
# O aumento da espessura (e) introduz uma nova camada de material com baixa 
# condutividade térmica (k). Esta camada atua como uma barreira física que 
# aumenta drasticamente a resistência térmica de condução (R_cond)[cite: 30]. Como a 
# perda de calor (q') é inversamente proporcional à resistência térmica total, o 
# aumento da espessura reduz a taxa de transferência de calor da tubulação para 
# o ambiente[cite: 30].

# b) Por que a minimização de CE(e), isoladamente, não representa o problema econômico completo?
# Se minimizássemos apenas o custo da energia perdida (CE), a matemática apontaria 
# para uma espessura de isolante tendendo ao infinito (e -> ∞), pois quanto mais 
# espesso for o isolante, menor será o calor perdido[cite: 30]. No entanto, isolantes 
# custam dinheiro[cite: 30]. O problema económico completo (C(e)) exige o equilíbrio 
# (trade-off) entre poupar energia e o custo de investimento (CI) de instalar e 
# comprar mais volume de material isolante[cite: 30].

# c) Qual a principal diferença entre a atualização do intervalo na razão áurea e a construção de novas aproximações pela interpolação quadrática?
# O método da Razão Áurea baseia-se numa redução geométrica constante (fator τ) do 
# intervalo de incerteza que contém o mínimo, avaliando apenas o sinal da função[cite: 31]. 
# Ele é lento (linear), mas incondicionalmente seguro e garantido a convergir[cite: 32]. A 
# Interpolação Quadrática, por outro lado, ajusta uma parábola usando os três 
# pontos atuais e salta diretamente para o vértice dessa parábola[cite: 32]. Esse salto pode 
# levar rapidamente ao mínimo real (convergência superlinear), mas é um método 
# arriscado que pode divergir ou exigir proteções computacionais se a curva não 
# tiver formato de parábola (A <= 0)[cite: 32].

# d) Qual método exigiu menos avaliações de C(e)? A resposta se manteve ao variar os pontos iniciais da interpolação?
# A Interpolação Quadrática geralmente exige muito menos avaliações da função C(e) 
# porque o salto direto para o vértice da parábola acelera a convergência[cite: 32]. 
# Contudo, esta eficiência extrema depende fortemente da escolha dos três pontos 
# iniciais (x1, x2, x3)[cite: 32]. Se os pontos iniciais estiverem muito longe do mínimo ou 
# não englobarem corretamente o "vale" da curva, o método quadrático pode acionar as 
# rotinas de proteção (fazendo passos de razão áurea) ou perder a eficiência. 
# A razão áurea tem um número de avaliações previsível e constante.

# e) Por que é importante impor limites físicos para a espessura do isolamento?
# Matematicamente, a função custo pode apresentar comportamentos indesejados para 
# valores extremos, como espessuras negativas (que não têm sentido físico)[cite: 30] ou 
# espessuras excessivamente grandes. Limites físicos (como [0, 0.12] m) garantem 
# que o algoritmo procure o mínimo numa zona onde a construção e instalação do 
# isolamento são tecnicamente viáveis para a tubulação[cite: 30], impedindo que as raízes 
# divirjam para soluções abstratas.

# f) Se a tubulação operasse por um tempo muito menor, como você esperaria que a espessura ótima se alterasse?
# Se o tempo de operação acumulado (t) fosse muito menor, a parcela relativa ao 
# custo da energia térmica perdida (CE) seria substancialmente menor no cômputo do 
# custo total (C)[cite: 30]. Sendo a energia gasta menor, o investimento financeiro para 
# comprar grandes espessuras de isolante deixaria de ser justificável. O mínimo 
# económico deslocar-se-ia, resultando numa espessura ótima menor (menos isolante)[cite: 30].

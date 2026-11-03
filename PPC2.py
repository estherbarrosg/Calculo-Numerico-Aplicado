import math
import cmath
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

def bairstow(a, r, s, tol=1e-5, max_iter=100):
    """
    Implementação do Método de Bairstow para encontrar raízes de um polinômio.
    a: lista de coeficientes em ordem decrescente [a_n, a_{n-1}, ..., a_0]
    r, s: chutes iniciais
    tol: tolerância para o erro relativo
    max_iter: número máximo de iterações
    """
    n = len(a) - 1
    raizes = []
    
    # Cópia dos coeficientes para podermos iterar (deflação)
    coef = a[:]
    
    while n >= 3:
        iteracao = 0
        erro_r = 1.0
        erro_s = 1.0
        
        while (erro_r > tol or erro_s > tol) and iteracao < max_iter:
            # Vetores b e c iniciados com zero
            b = [0.0] * (n + 1)
            c = [0.0] * (n + 1)
            
            # Cálculo dos b's (Divisão Sintética)
            b[0] = coef[0]
            b[1] = coef[1] + r * b[0]
            for i in range(2, n + 1):
                b[i] = coef[i] + r * b[i-1] + s * b[i-2]
                
            # Cálculo dos c's (Derivadas Parciais)
            c[0] = b[0]
            c[1] = b[1] + r * c[0]
            for i in range(2, n):
                c[i] = b[i] + r * c[i-1] + s * c[i-2]
                
            # Montagem do Sistema Linear 2x2 (Regra de Cramer)
            # c1*dr + c2*ds = -b0  --> (c[n-1]*dr + c[n-2]*ds = -b[n])
            # c2*dr + c3*ds = -b1  --> (c[n-2]*dr + c[n-3]*ds = -b[n-1])
            det_principal = c[n-1] * c[n-3] - c[n-2]**2
            
            if det_principal == 0:
                print("Determinante nulo. Mudando chute inicial...")
                r += 0.5
                s += 0.5
                continue
                
            dr = (-b[n] * c[n-3] - (-b[n-1]) * c[n-2]) / det_principal
            ds = (c[n-1] * (-b[n-1]) - c[n-2] * (-b[n])) / det_principal
            
            # Atualização dos valores
            r += dr
            s += ds
            
            # Cálculo dos erros relativos
            erro_r = abs(dr / r) if r != 0 else abs(dr)
            erro_s = abs(ds / s) if s != 0 else abs(ds)
            iteracao += 1
            
        # Extração do par de raízes usando a fórmula quadrática
        delta = r**2 + 4*s
        if delta >= 0:
            x1 = (r + math.sqrt(delta)) / 2.0
            x2 = (r - math.sqrt(delta)) / 2.0
        else:
            x1 = (r + cmath.sqrt(delta)) / 2.0
            x2 = (r - cmath.sqrt(delta)) / 2.0
            
        raizes.extend([x1, x2])
        
        # Deflação polinomial: O novo polinômio são os coeficientes b (menos o resto)
        coef = b[:-2]
        n -= 2

    # Tratamento para polinômios restantes de grau 1 ou 2
    if n == 2:
        r = -coef[1] / coef[0]
        s = -coef[2] / coef[0]
        delta = r**2 + 4*s
        if delta >= 0:
            x1 = (r + math.sqrt(delta)) / 2.0
            x2 = (r - math.sqrt(delta)) / 2.0
        else:
            x1 = (r + cmath.sqrt(delta)) / 2.0
            x2 = (r - cmath.sqrt(delta)) / 2.0
        raizes.extend([x1, x2])
    elif n == 1:
        x1 = -coef[1] / coef[0]
        raizes.append(x1)
        
    return raizes

# =====================================================================
# BLOCO DE EXECUÇÃO: Resolução do Sistema do APC2
# =====================================================================
if __name__ == "__main__":
    # O polinômio do APC2: P(λ) = λ^4 + 5λ^3 + 28λ^2 + 30λ + 75 = 0
    # Coeficientes em ordem decrescente de grau (do λ^4 para o termo independente)
    coef_apc2 = [1, 5, 28, 30, 75]
    
    # Chutes iniciais escolhidos (podem ser alterados para testar convergência)
    r_inicial = 1.0
    s_inicial = 1.0
    
    print("Iniciando Método de Bairstow para o sistema físico (APC2)...")
    raizes_sistema = bairstow(coef_apc2, r_inicial, s_inicial)
    
    print("\nAs raízes (autovalores) do sistema são:")
    for i, raiz in enumerate(raizes_sistema):
        # Formatando a saída para ficar elegante no terminal
        if isinstance(raiz, complex):
            print(f"λ{i+1} = {raiz.real:.4f} {'+' if raiz.imag >= 0 else '-'} {abs(raiz.imag):.4f}i")
        else:
            print(f"λ{i+1} = {raiz:.4f}")
            
            
            
        # =====================================================================
# ETAPA 2: Validação com Polinômio de 7ª Ordem
# =====================================================================
# Polinômio construído a partir das raízes: 0, 1, -1, 2, -2, 3, -3
# P(x) = x^7 - 14x^5 + 49x^3 - 36x = 0
coef_7 = [1, 0, -14, 0, 49, 0, -36, 0]

print("\n" + "="*50)
print("Iniciando validação com Polinômio de 7ª Ordem...")
raizes_7 = bairstow(coef_7, 1.0, 1.0)

print("Raízes encontradas (devem convergir para 0, 1, -1, 2, -2, 3, -3):")
for i, r in enumerate(raizes_7):
    # Se a parte imaginária for muito próxima de zero, imprimimos apenas a parte real
    if isinstance(r, complex) and abs(r.imag) > 1e-5:
        print(f"λ{i+1} = {r.real:.4f} {'+' if r.imag >= 0 else '-'} {abs(r.imag):.4f}i")
    else:
        print(f"λ{i+1} = {r.real if isinstance(r, complex) else r:.4f}")

# =====================================================================
# ETAPAS 6 e 7: Varredura Bidimensional e Fractal de Bairstow
# =====================================================================
print("\n" + "="*50)
print("Gerando o Fractal de Bairstow para o sistema físico...")
print("Isso pode demorar alguns segundos, aguarde a janela do gráfico abrir!")

def iteracoes_bairstow(a, r, s, tol=1e-5, max_iter=50):
    """Função enxuta para uso exclusivo na geração do Fractal (retorna iterações)"""
    n = len(a) - 1
    coef = a[:]
    iteracoes_totais = 0
    
    while n >= 3:
        iteracao = 0
        erro_r, erro_s = 1.0, 1.0
        
        while (erro_r > tol or erro_s > tol) and iteracao < max_iter:
            b = [0.0] * (n + 1)
            c = [0.0] * (n + 1)
            
            b[0] = coef[0]
            b[1] = coef[1] + r * b[0]
            for i in range(2, n + 1):
                b[i] = coef[i] + r * b[i-1] + s * b[i-2]
                
            c[0] = b[0]
            c[1] = b[1] + r * c[0]
            for i in range(2, n):
                c[i] = b[i] + r * c[i-1] + s * c[i-2]
                
            det_principal = c[n-1] * c[n-3] - c[n-2]**2
            if det_principal == 0:
                return max_iter # Falha na convergência retorna teto de iterações
                
            dr = (-b[n] * c[n-3] + b[n-1] * c[n-2]) / det_principal
            ds = (-c[n-1] * b[n-1] + c[n-2] * b[n]) / det_principal
            
            r += dr
            s += ds
            erro_r = abs(dr / r) if r != 0 else abs(dr)
            erro_s = abs(ds / s) if s != 0 else abs(ds)
            iteracao += 1
            
        iteracoes_totais += iteracao
        coef = b[:-2]
        n -= 2
        
    return iteracoes_totais

# Malha de varredura bidimensional no plano (r0, s0)
resolucao = 200 # Se o cálculo travar o seu computador, mude este valor para 100
r_vals = np.linspace(-3, 3, resolucao)
s_vals = np.linspace(-3, 3, resolucao)
fractal = np.zeros((resolucao, resolucao))

for i in range(resolucao):
    for j in range(resolucao):
        # Utiliza o polinômio original do APC2 para montar a matriz de convergência
        fractal[i, j] = iteracoes_bairstow([1, 5, 28, 30, 75], r_vals[j], s_vals[i])

# Plotagem do Fractal com mapa de cores
plt.figure(figsize=(8, 6))
plt.pcolormesh(r_vals, s_vals, fractal, cmap='turbo', shading='auto')
plt.colorbar(label='Número de Iterações até Convergir')
plt.title('Fractal de Bairstow - Zonas de Convergência (APC2)')
plt.xlabel('Chute Inicial (r_0)')
plt.ylabel('Chute Inicial (s_0)')
plt.tight_layout()
plt.show()
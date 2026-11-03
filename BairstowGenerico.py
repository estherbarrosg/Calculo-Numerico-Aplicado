import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PARTE 1: ALGORITMO GENÉRICO DO MÉTODO DE BAIRSTOW
# ==============================================================================
def bairstow_fator_quadratico(a, r0, s0, tol=1e-5, max_iter=50):
    """
    Executa o Método de Bairstow para extrair UM fator quadrático (x^2 - rx - s)
    do polinómio.
    
    Parâmetros:
    a        : Lista de coeficientes do polinómio (ordem crescente de grau).
               Ex: P(x) = a[0] + a[1]x + a[2]x^2 + ... + a[n]x^n
    r0, s0   : Estimativas iniciais para r e s.
    tol      : Tolerância de erro para o critério de convergência.
    max_iter : Limite de iterações.
    
    Retorna:
    (r, s, iteracoes, convergiu)
    """
    n = len(a) - 1  # Grau do polinómio
    
    if n < 3:
        raise ValueError("O polinómio deve ter no mínimo grau 3 para aplicar Bairstow desta forma.")
        
    r, s = r0, s0
    
    # Criamos os vetores 'b' (para o resto da divisão) e 'c' (para as derivadas parciais).
    # O tamanho é garantido ser no mínimo 4 para não dar erro de índice em c[3].
    tamanho_vetor = max(n + 1, 4)
    b = [0.0] * tamanho_vetor
    c = [0.0] * tamanho_vetor
    
    for iteracao in range(max_iter):
        
        # ----------------------------------------------------------------------
        # ETAPA A: Divisão Sintética Principal (Cálculo dos b_i)
        # Consiste em dividir P(x) por (x^2 - rx - s).
        # ----------------------------------------------------------------------
        b[n] = a[n]
        b[n-1] = a[n-1] + r * b[n]
        
        # O laço desce até ao termo independente (índice 0)
        for i in range(n-2, -1, -1):
            b[i] = a[i] + r * b[i+1] + s * b[i+2]
            
        # O resto da divisão está contido em b[1] e b[0].
        # O objetivo do método é fazer b[1] = 0 e b[0] = 0.
        
        # ----------------------------------------------------------------------
        # ETAPA B: Segunda Divisão Sintética (Cálculo das derivadas c_i)
        # Consiste em dividir o quociente anterior pelo mesmo fator quadrático.
        # ----------------------------------------------------------------------
        c[n] = b[n]
        c[n-1] = b[n-1] + r * c[n]
        
        # Desce apenas até c[1], pois o sistema Jacobiano só precisa de c[1], c[2] e c[3].
        for i in range(n-2, 0, -1):
            c[i] = b[i] + r * c[i+1] + s * c[i+2]
            
        # ----------------------------------------------------------------------
        # ETAPA C: Montagem e Resolução do Sistema Jacobiano Linear 2x2
        # O sistema é definido por:
        # [ c[2]  c[3] ] [ dr ]   [ b[1] ]
        # [ c[1]  c[2] ] [ ds ] = [ b[0] ]
        # ----------------------------------------------------------------------
        det = c[2]*c[2] - c[1]*c[3]
        
        if det == 0:
            # Se o determinante for zero, a matriz Jacobiana é singular (falha do método)
            return r, s, max_iter, False
            
        # Resolução utilizando a Regra de Cramer
        dr = (b[1]*c[2] - b[0]*c[3]) / det
        ds = (b[0]*c[2] - b[1]*c[1]) / det
        
        # ----------------------------------------------------------------------
        # ETAPA D: Atualização das estimativas
        # ----------------------------------------------------------------------
        r += dr
        s += ds
        
        # ----------------------------------------------------------------------
        # ETAPA E: Critério de Convergência
        # ----------------------------------------------------------------------
        erro_r = abs(dr / r) if r != 0 else abs(dr)
        erro_s = abs(ds / s) if s != 0 else abs(ds)
        
        if erro_r < tol and erro_s < tol:
            return r, s, iteracao + 1, True  # Convergiu!
            
    # Se sair do laço sem retornar, não convergiu no limite de iterações
    return r, s, max_iter, False


# ==============================================================================
# PARTE 2: VARREDURA BIDIMENSIONAL (GERAÇÃO DO FRACTAL DE BAIRSTOW)
# ==============================================================================
def varredura_bairstow_fractal(a, r_bounds, s_bounds, resolucao=300, max_iter=40, tol=1e-4):
    """
    Mapeia a convergência do método num plano bidimensional de estimativas (r0, s0).
    
    r_bounds : Tupla (r_min, r_max) delimitando o eixo x.
    s_bounds : Tupla (s_min, s_max) delimitando o eixo y.
    resolucao: Quantidade de pontos em cada eixo (Ex: 300 gera uma malha 300x300).
    """
    print("A iniciar varredura sistemática... (isto pode demorar alguns segundos dependendo da resolução)")
    
    # 1. Criação dos vetores de varredura uniformemente espaçados
    r_vals = np.linspace(r_bounds[0], r_bounds[1], resolucao)
    s_vals = np.linspace(s_bounds[0], s_bounds[1], resolucao)
    
    # Matriz para armazenar o número de iterações (responsável pelas cores do fractal)
    mapa_fractal = np.zeros((resolucao, resolucao))
    
    # 2. Laço duplo varrendo todo o plano (r, s)
    for i, r0 in enumerate(r_vals):
        for j, s0 in enumerate(s_vals):
            
            # Executa Bairstow para o par (r0, s0) atual
            _, _, iteracoes, convergiu = bairstow_fator_quadratico(
                a, r0, s0, tol=tol, max_iter=max_iter
            )
            
            # Se não convergir, pintamos com a cor de "fundo" mapeando-o para max_iter
            if not convergiu:
                mapa_fractal[j, i] = max_iter
            else:
                # O índice [j, i] obedece à lógica matricial: j (linhas) = eixo s, i (colunas) = eixo r
                mapa_fractal[j, i] = iteracoes
                
    # 3. Plotagem Visual do Fractal usando Matplotlib
    plt.figure(figsize=(10, 8))
    
    # O cmap='inferno' ou 'magma' evidencia muito bem a profundidade das iterações
    plt.imshow(mapa_fractal, extent=[r_bounds[0], r_bounds[1], s_bounds[0], s_bounds[1]], 
               origin='lower', cmap='inferno', aspect='auto')
               
    barra_cores = plt.colorbar(label='Nº de Iterações para Convergência')
    
    # Formatações do gráfico
    plt.title('Fractal de Bairstow: Regiões de Convergência', fontsize=14, fontweight='bold')
    plt.xlabel('Estimativa Inicial $r_0$', fontsize=12)
    plt.ylabel('Estimativa Inicial $s_0$', fontsize=12)
    plt.grid(False)
    plt.tight_layout()
    
    # Salva o resultado e exibe no ecrã (perfeito para colocares no teu repositório)
    plt.savefig('fractal_bairstow.png', dpi=300)
    print("Gráfico 'fractal_bairstow.png' gerado com sucesso!")
    plt.show()

# ==============================================================================
# BLOCO PRINCIPAL DE TESTE (Exemplo Prático)
# ==============================================================================
if __name__ == "__main__":
    
    # Exemplo: Polinómio P(x) = x^5 - 1
    # Raízes são as raízes complexas da unidade. 
    # Ordem dos coeficientes: a[0] + a[1]x + ... + a[5]x^5
    coeficientes = [-1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    
    # Define a janela do plano (r0, s0) para a varredura
    r_min, r_max = -2.0, 2.0
    s_min, s_max = -2.0, 2.0
    
    # Executa a geração do fractal (Resolução 300x300 = 90.000 sistemas resolvidos)
    # Podes baixar a resolução para 100 se quiseres que execute quase instantaneamente.
    varredura_bairstow_fractal(coeficientes, (r_min, r_max), (s_min, s_max), resolucao=300)
# PPC 01 - Método de Runge-Kutta na Sedimentação de Partículas

### Resumo Operacional
Implementação computacional do método de Runge-Kutta de quarta ordem (RK4) para resolver a equação diferencial ordinária que descreve o processo de sedimentação de uma esfera num fluido viscoso. O código avalia os efeitos inerciais do escoamento modelados pela correção de Oseen.

### Dicionário de Variáveis
* `St`: Número de Stokes (Adimensional, float) - Mede a escala de tempo de relaxação.
* `Re_s`: Número de Reynolds da partícula (Adimensional, float).
* `h`: Passo de tempo de integração (Adimensional, float).
* `t_atual` e `v_atual`: Variáveis de estado do tempo e da velocidade (Adimensional, float).
* `lista_t` e `lista_v`: Estruturas para armazenamento dinâmico do histórico (array).

### Dependências e Bibliotecas
* `numpy`: Utilizado estritamente para funcionalidades matemáticas básicas (manipulação de matrizes e vetores).
* `matplotlib.pyplot`: Utilizado exclusivamente para a geração de gráficos analíticos.

### Especificações de I/O (Entradas e Saídas)
* **Inputs:** Os parâmetros de modelagem estão configurados nas funções de cada análise no corpo principal do script.
* **Outputs:** Gráficos salvos localmente em `.png` dentro do diretório `outputs`, plotando o histórico temporal da velocidade $v(t)$.

### Procedimentos de Execução
Para executar o script via interpretador Python, abra o terminal na pasta raiz e utilize o comando:
```bash
python ppc1_esther_gomes.py
```

### Validação Metodológica
O algoritmo foi submetido a uma análise comparativa demonstrando sua precisão. A implementação do RK4 recuperou perfeitamente as soluções analíticas exatas de referência para o caso limite sem inércia ($Re \to 0$, por Fatores Integrantes) e para o modelo não linear com arrasto de Oseen ($Re \neq 0$, Equação de Riccati).

### Bibliografia Específica
1. SOBRAL, Y. D.; OLIVEIRA, T. F.; CUNHA, F. R. "On the unsteady forces during the motion of a sedimenting particles". Powder Technology, 178 (2007), 129-141.

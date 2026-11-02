# PPC 01 - Métodos Intervalares para Raízes de Funções (Bisseção e Falsa Posição)

> **Status:** método da **bisseção** implementado e validado. O método da **falsa posição** será adicionado nesta mesma pasta (`falsa_posicao.py`).

## 1. Resumo operacional

Escoamento uniforme de água em um canal aberto retangular de largura $B$, inclinado em relação à horizontal. Dada a vazão volumétrica $Q$, o programa determina a **profundidade $H$ da lâmina d'água** e, a partir dela, a **velocidade média $U$** do escoamento.

A relação entre $Q$ e $H$ é não linear (equação de Manning), então $H$ não pode ser isolado algebricamente. O problema é reescrito como a busca do zero de uma função, $f(H) = 0$, e resolvido numericamente pelo método da bisseção.

Parâmetros do problema:

| Parâmetro | Símbolo | Valor | Unidade |
|---|---|---|---|
| Declividade do canal | $S$ | 0,0002 | m/m |
| Largura do canal | $B$ | 20 | m |
| Coeficiente de rugosidade de Manning | $n$ | 0,03 | s/m^(1/3) |
| Vazão prescrita | $Q$ | 5 | m³/s |

## 2. Fundamentação matemática

A vazão que atravessa a seção do canal é

$$Q = U A_c = U B H$$

com $A_c = BH$ a área molhada. A velocidade média é estimada pela equação de Manning:

$$U = \frac{1}{n} R^{2/3} S^{1/2}, \qquad R = \frac{A_c}{P} = \frac{BH}{B + 2H}$$

em que $R$ é o raio hidráulico e $P = B + 2H$ é o perímetro molhado. Combinando as expressões:

$$Q(H) = \frac{\sqrt{S}}{n} \, \frac{(BH)^{5/3}}{(B + 2H)^{2/3}}$$

A profundidade procurada é a raiz de

$$f(H) = Q(H) - Q = 0$$

Fisicamente, a raiz é a profundidade para a qual a vazão calculada por Manning é igual à vazão prescrita.

## 3. Algoritmo implementado (bisseção)

O método parte de um intervalo $[H_l, H_u]$ em que $f$ muda de sinal, isto é, $f(H_l)\,f(H_u) < 0$. Para uma função contínua, o Teorema de Bolzano garante pelo menos uma raiz no intervalo.

1. Intervalo inicial: $H_l = 0{,}5$ m e $H_u = 1{,}0$ m. Verifica-se $f(H_l)\,f(H_u) < 0$.
2. Em cada iteração $k$, calcula-se o ponto médio $H_r = \dfrac{H_l + H_u}{2}$ e o valor $f(H_r)$.
3. **Critério de parada:** $\lvert f(H_r) \rvert < 10^{-6}$ (limitado a 50 iterações).
4. Atualização do intervalo, mantendo sempre a mudança de sinal:
   - se $f(H_l)\,f(H_r) < 0$, então $H_u = H_r$;
   - caso contrário, $H_l = H_r$.
5. Após a convergência, calcula-se $U = \dfrac{Q}{BH}$.

A tabela também mostra o erro relativo aproximado entre iterações consecutivas, apenas como informação sobre a convergência (não é usado como critério de parada):

$$\varepsilon_a = \left| \frac{H_r^{\text{novo}} - H_r^{\text{anterior}}}{H_r^{\text{novo}}} \right| \times 100\%$$

## 4. Organização da pasta

```
PPC_01/
├── README.md              Este arquivo
├── problema_canal.py      Parâmetros e equações do problema físico (Manning); sem método numérico
├── bissecao.py            Método da bisseção (função genérica) e execução para o problema do canal
└── outputs/
    └── resultados_bissecao.txt   Relatório gerado pelo programa (tabela de iterações e resultados)
```

`problema_canal.py` é separado de propósito: o mesmo problema será resolvido pela falsa posição importando as mesmas funções, sem repetir código.

## 5. Dicionário de variáveis

### `problema_canal.py`

| Variável no código | Significado | Unidade | Tipo |
|---|---|---|---|
| `declividade` | $S$, inclinação do canal | m/m | float |
| `largura_canal` | $B$, largura do canal | m | float |
| `rugosidade_manning` | $n$, coeficiente de rugosidade | s/m^(1/3) | float |
| `vazao_conhecida` | $Q$, vazão prescrita | m³/s | float |
| `profundidade` | $H$, profundidade da lâmina d'água | m | float |
| `area_molhada` | $A_c = BH$, área da seção ocupada pela água | m² | float |
| `perimetro_molhado` | $P = B + 2H$ | m | float |
| `coeficiente_manning` | $\sqrt{S}/n$ | m^(1/3)/s | float |

### `bissecao.py`

| Variável no código | Significado | Unidade | Tipo |
|---|---|---|---|
| `limite_inferior_inicial` | $H_l$ inicial | m | float |
| `limite_superior_inicial` | $H_u$ inicial | m | float |
| `tolerancia` | tolerância do critério de parada, $10^{-6}$ | m³/s | float |
| `max_iteracoes` | número máximo de iterações | - | int |
| `limite_inferior`, `limite_superior` | $H_l$ e $H_u$ atuais, atualizados a cada iteração | m | float |
| `f_inferior`, `f_superior` | $f(H_l)$ e $f(H_u)$ | m³/s | float |
| `aproximacao` | $H_r$, ponto médio do intervalo | m | float |
| `f_aproximacao` | $f(H_r)$, resíduo da equação | m³/s | float |
| `aproximacao_anterior` | $H_r$ da iteração anterior | m | float ou None |
| `erro_relativo_percentual` | $\varepsilon_a$ entre iterações consecutivas | % | float ou None |
| `iteracao` | contador $k$ | - | int |
| `historico` | uma tupla por iteração: (k, $H_l$, $H_u$, $H_r$, $f(H_r)$, $\varepsilon_a$) | - | list de tuplas |
| `convergiu` | indica se o critério de parada foi atendido | - | bool |
| `profundidade` (em `main`) | raiz encontrada, $H$ | m | float ou None |
| `velocidade` | $U$, velocidade média | m/s | float |
| `relatorio` | linhas de texto do relatório | - | list de str |

## 6. Dependências e bibliotecas

Apenas a biblioteca padrão do Python 3 (`math` e `os`). Não é necessário instalar nenhum pacote, e nenhuma função pronta de busca de raízes é usada.

## 7. Entradas e saídas (I/O)

**Entradas:** não há arquivos de entrada. Os parâmetros do problema ficam no início de `problema_canal.py` e a configuração do método (intervalo, tolerância, máximo de iterações) fica no início de `bissecao.py`.

**Saídas:**
- Console: tabela com $k$, $H_l$, $H_u$, $H_r$, $f(H_r)$ e $\varepsilon_a$ por iteração, seguida da profundidade encontrada, número de iterações, verificação da vazão e velocidade média.
- Arquivo `outputs/resultados_bissecao.txt`: o mesmo relatório do console.

## 8. Procedimentos de execução

Dentro da pasta `PPC_01`:

```bash
python bissecao.py
```

Ou a partir da raiz do repositório (o arquivo de saída é salvo sempre em `PPC_01/outputs/`):

```bash
python PPC_01/bissecao.py
```

## 9. Resultados

Intervalo inicial $[0{,}5;\ 1{,}0]$ m, com $f(0{,}5) = -2{,}125377$ e $f(1{,}0) = 3{,}847664$ (mudança de sinal confirmada). Relatório completo em [`outputs/resultados_bissecao.txt`](outputs/resultados_bissecao.txt).

| Grandeza | Valor |
|---|---|
| Profundidade $H$ | 0,702293 m |
| Iterações até a convergência | 22 |
| Velocidade média $U$ | 0,3560 m/s |

## 10. Validação metodológica

**a) Resíduo da equação de Manning.** Substituindo $H = 0{,}702293$ m de volta na equação de Manning, o programa obtém $Q(H) = 5{,}000000$ m³/s, igual à vazão prescrita. Na última iteração, $f(H_r) = 2{,}37 \times 10^{-7}$ m³/s, abaixo da tolerância de $10^{-6}$.

**b) Comportamento teórico da bisseção.** A largura do intervalo cai pela metade a cada iteração. Na iteração 22, a largura é $\dfrac{0{,}5}{2^{21}} \approx 2{,}4 \times 10^{-7}$ m, o que garante que o erro de $H_r$ em relação à raiz exata é menor que cerca de $1{,}2 \times 10^{-7}$ m. As colunas $H_l$ e $H_u$ da tabela mostram esse comportamento.

**c) Comparação com solução analítica aproximada.** Para um canal muito largo ($B \gg H$), vale $P \approx B$, e a equação de Manning pode ser invertida em forma fechada:

$$H \approx \left( \frac{Q\,n}{B\sqrt{S}} \right)^{3/5} \approx 0{,}6835 \text{ m}$$

O valor numérico, 0,702293 m, é cerca de 2,7% maior. A diferença é coerente com a aproximação: ao desprezar o termo $2H$ no perímetro molhado, o raio hidráulico fica superestimado e a profundidade necessária para transportar a mesma vazão fica subestimada.

**d) Comparação com a falsa posição.** *A ser preenchida após a implementação do segundo método.*

## 11. Bibliografia específica

- CHAPRA, Steven C.; CANALE, Raymond P. *Métodos numéricos para engenharia*. McGraw-Hill. (Capítulo sobre métodos intervalares: bisseção e falsa posição.)
- UNIVERSIDADE DE BRASÍLIA. Departamento de Engenharia Mecânica. *Atividade Prática de Programação: Métodos da Bisseção e da Falsa Posição*. Cálculo Numérico Aplicado (ENM0227), 2026.

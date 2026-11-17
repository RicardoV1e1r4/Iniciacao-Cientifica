# Série de Taylor de ordem qualquer para a amplitude de um caminho

Este documento explica as funções `taylorAmplitude()` e `chebyshevPoints()` do arquivo `underwaterFunctions.py`. Juntas, elas aproximam a variação da amplitude de **um** caminho do canal acústico submarino dentro de um bloco de duração $T_{bl}$, usando uma série de Taylor de ordem $N$ qualquer.

O programa `testeFuncoes07.py` mostra as duas funções em uso.

## Sumário

1. [Motivação e base teórica](#1-motivação-e-base-teórica)
2. [Função `chebyshevPoints()`](#2-função-chebyshevpoints)
3. [Função `taylorAmplitude()`](#3-função-tayloramplitude)
4. [Exemplo de uso](#4-exemplo-de-uso)
5. [Validação](#5-validação)
6. [Limitações e cuidados](#6-limitações-e-cuidados)
7. [Referências](#7-referências)

---

## 1. Motivação e base teórica

O livro de Zhou e Wang (2014) modela o canal como uma soma de $N_{pa}$ caminhos discretos (eq. 1.11):

$$
h(t;\tau) = \sum_{p=1}^{N_{pa}} A_p(t)\,\delta\big(\tau - \tau_p(t)\big)
$$

em que $A_p(t)$ é a amplitude e $\tau_p(t)$ é o atraso do caminho $p$.

A hipótese mais simples considera a amplitude constante dentro do bloco (eq. 1.12). A seção 1.3.4 generaliza essa hipótese e aproxima a amplitude por um polinômio de Taylor de ordem $N_{amp}$ (eq. 1.19):

$$
A_p(t) \approx \sum_{n=0}^{N_{amp}} \frac{(-1)^n}{n!}\,A_p^{(n)}\,t^n
$$

As funções deste documento implementam essa ideia para qualquer ordem. Na prática, elas precisam responder a uma pergunta: **como obter as derivadas $A_p^{(n)}$ a partir de valores da amplitude medidos dentro do bloco?**

A resposta usada aqui é a interpolação polinomial. Um polinômio de grau $N$ fica totalmente determinado por $N+1$ pontos. Então, medindo o ganho em $N+1$ instantes do bloco, é possível encontrar o único polinômio de grau $N$ que passa por todos eles e, a partir dos coeficientes, calcular as derivadas.

> **Por que a ordem N exige N + 1 pontos?** Uma reta (ordem 1) precisa de 2 pontos. Uma parábola (ordem 2) precisa de 3. Em geral, um polinômio de ordem $N$ tem $N+1$ coeficientes desconhecidos, e cada ponto medido fornece uma equação.

---

## 2. Função `chebyshevPoints()`

### Para que serve

Ela escolhe **onde** amostrar a amplitude dentro do bloco. A escolha dos instantes muda muito a qualidade da aproximação quando a ordem é alta.

### Assinatura

```python
uwf.chebyshevPoints(M, Tbl, t_start=0.0, lobatto=False)
```

| Parâmetro | Significado |
|---|---|
| `M` | Número de pontos. Para uma série de ordem $N$, use `M = N + 1`. |
| `Tbl` | Duração do bloco, em segundos. |
| `t_start` | Instante inicial do bloco, em segundos. O padrão é 0. |
| `lobatto` | Escolhe a variante dos pontos. Veja abaixo. |

**Retorno:** um vetor com os `M` instantes, em ordem crescente, dentro do intervalo $[t_{start},\ t_{start} + T_{bl}]$.

### Como funciona

**Passo 1. Pontos no intervalo [−1, 1].** A função calcula os pontos de Chebyshev no intervalo padrão.

Com `lobatto=False`, são os pontos de Chebyshev clássicos:

$$
x_k = \cos\left(\frac{(2k-1)\,\pi}{2M}\right), \qquad k = 1, 2, \dots, M
$$

Com `lobatto=True`, são os pontos de Chebyshev–Lobatto:

$$
x_k = \cos\left(\frac{k\,\pi}{M-1}\right), \qquad k = 0, 1, \dots, M-1
$$

**Passo 2. Conversão para o bloco.** Cada ponto é levado para o intervalo de tempo do bloco:

$$
t_k = t_{start} + \frac{T_{bl}}{2}\,(x_k + 1)
$$

**Passo 3. Ordenação.** Os cossenos saem do maior para o menor, então o vetor é ordenado em ordem crescente de tempo.

### Interpretação geométrica

Marque pontos igualmente espaçados ao longo de um semicírculo e projete-os sobre o diâmetro. As projeções são os pontos de Chebyshev. Como o arco quase não avança na horizontal perto das pontas, os pontos se acumulam nas bordas do bloco.

### Por que usar esses pontos

O erro de um polinômio interpolador de grau $N$ é:

$$
f(t) - P_N(t) = \frac{f^{(N+1)}(\xi)}{(N+1)!}\,\prod_{k=0}^{N}(t - t_k)
$$

O primeiro fator depende da física do problema. O segundo, o produto, depende **apenas** de onde os pontos foram colocados.

- **Pontos igualmente espaçados.** O produto é pequeno no meio do bloco, mas cresce muito perto das bordas. Em ordens altas, isso causa oscilações perto das bordas, o chamado **fenômeno de Runge**.
- **Pontos de Chebyshev.** Eles são justamente os que tornam o maior valor desse produto o menor possível no intervalo inteiro. O erro fica distribuído de forma equilibrada.

### Clássicos ou Lobatto?

| Variante | Inclui as bordas do bloco? | Quando usar |
|---|---|---|
| Clássicos (`lobatto=False`) | Não | Quando se quer o menor erro máximo possível. |
| Lobatto (`lobatto=True`) | Sim | Quando a curva precisa passar exatamente pelo início e pelo fim do bloco. |

Exemplo com $T_{bl} = 60$ s e `M = 5`:

| Tipo de ponto | Instantes (s) |
|---|---|
| Igualmente espaçados | 0; 15; 30; 45; 60 |
| Chebyshev clássicos | 1,47; 12,37; 30; 47,63; 58,53 |
| Chebyshev–Lobatto | 0; 8,79; 30; 51,21; 60 |

---

## 3. Função `taylorAmplitude()`

### Para que serve

Ela calcula a amplitude $A_p(t)$ de um caminho ao longo de todo o bloco, a partir de alguns ganhos medidos, usando uma série de Taylor de ordem $N$.

### Assinatura

```python
uwf.taylorAmplitude(A0, gains, t_samples, time_grid, order=None, return_derivatives=False)
```

| Parâmetro | Significado |
|---|---|
| `A0` | Amplitude de referência do caminho, em escala linear. Em geral, a amplitude no centro do bloco. |
| `gains` | Ganhos relativos $g_k = A(t_k)/A_0$ medidos nos instantes `t_samples`. |
| `t_samples` | Instantes, em segundos, em que os ganhos foram medidos. Podem vir de `chebyshevPoints()` ou de qualquer outra escolha. |
| `time_grid` | Vetor de tempo onde a amplitude será calculada. O bloco vai de `time_grid[0]` a `time_grid[-1]`. |
| `order` | Ordem $N$ da série. Se for `None`, usa $N$ igual ao número de ganhos menos 1. |
| `return_derivatives` | Se `True`, também devolve as derivadas da amplitude no centro do bloco. |

**Retorno:**

- `A_t`: vetor com a amplitude $A_p(t)$ em cada instante de `time_grid`.
- `derivatives` (só com `return_derivatives=True`): vetor com $A^{(0)}, A^{(1)}, \dots, A^{(N)}$ no centro do bloco, em unidades de amplitude por $s^n$.

### Como funciona

**Passo 1. Verificação das entradas.** A função confere se `gains` e `t_samples` têm o mesmo tamanho e se há pelo menos $N+1$ ganhos. Caso contrário, gera um erro explicando o problema:

```
ValueError: A ordem 3 exige pelo menos 4 ganhos; foram dados 2.
```

**Passo 2. Normalização do tempo.** O tempo é centrado no meio do bloco e dividido pela meia-duração $h$:

$$
t_c = \frac{t_{inicio} + t_{fim}}{2}, \qquad h = \frac{t_{fim} - t_{inicio}}{2}, \qquad u = \frac{t - t_c}{h}
$$

Assim, $u = -1$ no início do bloco, $u = 0$ no centro e $u = +1$ no fim.

> **Por que normalizar?** Sem normalização, um bloco de 60 s levaria a termos como $30^8 \approx 6{,}6 \times 10^{11}$, enquanto os coeficientes correspondentes seriam minúsculos. Misturar números tão diferentes faz o computador perder precisão. Com $u$ entre −1 e 1, todas as potências $u^n$ também ficam entre −1 e 1.

**Passo 3. Ajuste do polinômio.** A função encontra os coeficientes $c_n$ do polinômio

$$
g(u) = \sum_{n=0}^{N} c_n\,u^n
$$

com `np.polynomial.polynomial.polyfit`.

- Com exatamente $N+1$ ganhos, o polinômio **passa por todos eles** (interpolação).
- Com mais de $N+1$ ganhos, o ajuste é feito por **mínimos quadrados**. Isso é útil quando os ganhos têm ruído, como em um receptor real.

**Passo 4. Cálculo da amplitude.** O polinômio é avaliado em todos os instantes do bloco e multiplicado pela amplitude de referência:

$$
A_p(t) = A_0 \sum_{n=0}^{N} c_n\,u(t)^n
$$

**Passo 5. Derivadas (opcional).** Pela regra da cadeia, como $du/dt = 1/h$, a derivada de ordem $n$ no centro do bloco é:

$$
A^{(n)} = \frac{d^n A_p}{dt^n}\bigg|_{t = t_c} = A_0\,\frac{n!\;c_n}{h^n}
$$

Com essas derivadas, a amplitude pode ser escrita exatamente na forma de série de Taylor em torno do centro do bloco:

$$
A_p(t) = \sum_{n=0}^{N} \frac{A^{(n)}}{n!}\,(t - t_c)^n
$$

### Relação com a equação 1.19 do livro

A ideia é a mesma, mas há duas diferenças de convenção:

| Aspecto | Livro (eq. 1.19) | `taylorAmplitude()` |
|---|---|---|
| Origem do tempo | Início do bloco, $t \in [0, T_{bl}]$ | Centro do bloco, $t - t_c \in [-T_{bl}/2,\ T_{bl}/2]$ |
| Sinal dos termos | Alternado, $(-1)^n$ | Sempre positivo |

As duas formas descrevem o mesmo tipo de polinômio. Só os valores dos coeficientes mudam. Expandir em torno do centro reduz o erro, porque o ponto mais distante da origem fica a apenas $T_{bl}/2$.

### Relação com a função antiga `taylorAmplitudeVariation()`

A função antiga, de ordem 1 ou 2, continua no arquivo. A nova é a generalização dela. Com ordem 2 e ganhos no início, no centro e no fim do bloco, as duas dão exatamente o mesmo resultado.

---

## 4. Exemplo de uso

Trecho baseado em `testeFuncoes07.py`:

```python
import numpy as np
import underwaterFunctions as uwf

Tbl   = 60                          # Duração do bloco em s
ordem = 4                           # Ordem da série de Taylor
M     = ordem + 1                   # Número de ganhos necessários
tempo = np.linspace(0, Tbl, 1001)   # Vetor de tempo do bloco

# 1. Escolhe os instantes de amostragem
t_k = uwf.chebyshevPoints(M, Tbl)

# 2. Mede a amplitude nesses instantes
#    (amplitude() é a função do teste que calcula a amplitude exata)
A0    = amplitude(Tbl / 2)          # Referência: centro do bloco
gains = amplitude(t_k) / A0         # Ganhos relativos

# 3. Calcula a amplitude em todo o bloco e as derivadas no centro
A_t, derivadas = uwf.taylorAmplitude(A0, gains, t_k, tempo,
                                     order=ordem, return_derivatives=True)
```

---

## 5. Validação

As funções foram verificadas de três formas:

1. **Polinômio conhecido.** Com ganhos tirados de um polinômio cúbico, a função de ordem 3 reproduziu a curva exatamente e devolveu as derivadas corretas.
2. **Compatibilidade.** Com ordem 2 e três pontos igualmente espaçados, o resultado coincidiu com o da função antiga `taylorAmplitudeVariation()`.
3. **Erros de uso.** Com menos ganhos do que a ordem exige, a função gera um erro com uma mensagem clara.

O `testeFuncoes07.py` compara a aproximação com a amplitude exata de um receptor que anda a 5 m/s em um bloco de 60 s. Erro relativo máximo em cada ordem:

| Ordem | Pontos de Chebyshev | Pontos igualmente espaçados |
|---|---|---|
| 1 | 1,03 % | 1,31 % |
| 2 | 0,27 % | 0,41 % |
| 4 | 0,0024 % | 0,0040 % |
| 6 | 0,00009 % | 0,00025 % |

Os pontos de Chebyshev dão erro menor em todas as ordens.

---

## 6. Limitações e cuidados

- **Variações bruscas.** A série funciona bem quando a amplitude varia suavemente dentro do bloco. Se o receptor passa muito perto do transmissor no meio do bloco, a amplitude forma um pico estreito e são necessárias ordens altas. Nesse caso, os pontos igualmente espaçados podem piorar com o aumento da ordem, por causa do fenômeno de Runge. Os pontos de Chebyshev continuam convergindo.
- **Pontos clássicos não incluem as bordas.** Com `lobatto=False`, a curva não passa exatamente pelo início e pelo fim do bloco. Use `lobatto=True` se isso for necessário.
- **Blocos reais são curtos.** Blocos OFDM reais duram entre 0,1 e 0,3 s. Nesse intervalo, a amplitude de grande escala quase não varia, e a ordem 1 ou 2 já basta. Os blocos de 60 s dos testes servem apenas para tornar a curvatura visível.
- **Estimação em um receptor real.** Nas simulações, os ganhos são calculados com exatidão. Em um receptor real, cada coeficiente precisa ser estimado a partir do sinal recebido, com ruído. Cada ordem a mais é um parâmetro a mais para estimar, e por isso o livro trabalha com ordens baixas.
- **Um caminho por vez.** As funções tratam um único caminho. Para o canal completo, é preciso aplicá-las a cada caminho, com os seus próprios ganhos, e somar as contribuições conforme a eq. 1.11.
- **Ordens muito altas.** A normalização do tempo mantém o cálculo estável em ordens moderadas. Em ordens muito altas, o ajuste pode perder precisão, e a matriz do ajuste fica mal condicionada.

---

## 7. Referências

- S. Zhou e Z. Wang, *OFDM for Underwater Acoustic Communications*. Wiley, 2014. Seção 1.3, equações 1.11 a 1.19.
- X. Xu, Z.-H. Wang, S. Zhou e L. Wan, "Parameterizing both path amplitude and delay variations of underwater acoustic channels for block decoding of orthogonal frequency division multiplexing", *J. Acoust. Soc. Am.*, pp. 4672–4679, 2012. Citado no livro acima como referência [442].

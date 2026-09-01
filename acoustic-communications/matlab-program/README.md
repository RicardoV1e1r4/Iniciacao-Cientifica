# Visão Geral do programa
O objetivo do código é implementar um **modelo estatístico de canal acústico submarino (UWA — Underwater Acoustic)**. Sendo a ideia principal, gerar uma **resposta ao impulso de um canal multipercurso**, considerando:

1. **Distribuição aleatória dos atrasos** dos caminhos;
2. **Ganhos aleatórios** de cada caminho;
3. **Efeito Doppler** sobre os atrasos;
4. Montagem de **resposta ao impulso discreta** $\large h$.

## 1. O que o programa está simulando?
Imagine que você queira transmitir um sinal acústico debaixo d'água:
```
             Transmissor
                  |
                  |\
                  | \ caminho 1
                  |  \________________ Receptor
                  |
                  |\
                  | \ caminho 2
                  |  \
                  |   \______________ Receptor
                  |
                  |\
                  | \ caminho 3
                  |  \
                  |   \______________ Receptor
```

Vê-se que é o mesmo receptor, mas o sinal transmitido pode ir por vários caminhos. O sinal sofre diversas reflexões. Entre elas, estão:
- fundo do mar;
- superfície;
- obstáculos;
- diferentes regiões da água.
Com isso, o receptor vai receber várias do mesmo sinal, só que, cada uma com:
- um _atraso diferente_;
- uma _amplitude diferente_;

O canal pode ser aproximado matematicamente pela expressão:
$$\large h(t) \space = \space \sum_{l=0}^{L-1}h_l \space \delta(t \space - \space\tau_l)$$
onde:
- $\large L =$ número de caminhos;
- $\large h_l =$ ganho do caminho $\large l$;
- $\large \tau_l =$ atraso do caminho $\large l$;
- $\large \delta(\cdot) =$ impulso.

O programa serve justamente para gerar valores para os parâmetros $\large h_l \text{ e } \tau_l$.

## 2. Entradas e saídas do programa
Declaração da função:
```Matlab
function [h,delay_bar,gain_tap,varargout] = acousticChannel(setup,gain,delay,doppler)
```
Parâmetros de entrada:
```
setup
gain
delay
doppler
```
Parâmetros de saída:
```
h
delay_bar
gain_tap
varargout
```

### 2.1 Entradas
#### setup
**setup.Ts** -> é o período de amostragem $\large T_s$;

**setup.paths** -> é o números de caminhos multipercurso;

**setup.delayspread** -> Representa o **espalhamento temporal do canal**.
Ele influencia a velocidade com que os ganhos dos caminhos diminuem à medida que o atraso aumenta.

#### gain
**gain.attenuation** -> é a atenuação do canal em dB.

#### delay
**delay.mean** -> Esse parâmetro é utilizado como média da distribuição exponencial dos atrasos.
Por exemplo:
```
delay.mean = 0.002;
```
significa uma média de: $\large 2ms$ para o intervalo entre caminhos.

#### doppler
**doppler.type** -> os tipos são: '_none_', '_uniform_' e '_non-uniform_'.
**doppler.velocity**

##### Sem usar o efeito Doppler
```
doppler.type = 'none';
```
##### Usando Doppler uniforme
```
doppler.type = 'uniform';
doppler.velocity = 10;
```
Nesse caso todos os caminhos sofrem o mesmo fator Doppler.
##### Usando Doppler não-uniforme
A intenção é que cada caminho tenha uma velocidade diferente:
```
doppler.type = 'non-uniform';
doppler.velocity = [5 8 10 12 ...];
```

## 3. Análise do código
### 3.1 Primeira etapa: geração dos atrasos
O código começa com:
```
distribution_delay = makedist('exponential','mu',delay.mean);
```
Aqui ele cria uma **distribuição exponencial**.
A distribuição exponencial é usada para modelar os intervalos entre os caminhos multipercurso.
Conceitualmente: $\Large \Delta \tau_l \space \thicksim \space Exponential(\mu)$
#### Gerando os intervalos aleatórios
Depois:
```
Dtau = random(distribution_delay,setup.paths,1);
```
Por exemplo, se:
```
setup.paths = 5;
```
podemos ter algo como:
```
Dtau =
    0.0008
    0.0015
    0.0003
    0.0021
    0.0010
```
Esses valores representam **intervalos entre caminhos**, não necessariamente os atrasos absolutos.

#### Convertendo os atrasos para tempo discreto
Depois:
```
Dtau_index = ceil(Dtau/setup.Ts);
```
Aqui acontece algo muito importante.
O código transforma o intervalo de tempo em número de amostras.
A equação é: $\Large \Delta n_l​ \space = \space \Big{[} \frac{\Delta \tau_l}{Ts} \Big{]}$
Por exemplo, se:
$\large Ts​=1ms$
e:
$\large \Delta \tau \space = \space 2.3ms$
então:
$\Large \Delta n​ \space = \space \Big{[} \frac{2.3}{1} \Big{]} \space = \space 2.3$
Ou seja, o caminho será colocado aproximadamente na amostra 3.

#### Calculando os atrasos acumulados
Depois:
```
delay_index = cumsum(Dtau_index);
```
`cumsum` faz uma soma acumulada.

Por exemplo:
```
Dtau_index =
    2
    1
    3
    2
```

vira:
```
delay_index =
    2
    3
    6
    8
```

Isso representa:
```
Caminho 1 → índice 2
Caminho 2 → índice 3
Caminho 3 → índice 6
Caminho 4 → índice 8
```

Portanto, os caminhos vão sendo posicionados ao longo do eixo temporal.

#### Convertendo novamente para tempo
Depois:
```
delay = (delay_index - 1)*setup.Ts;
```
Agora os índices são convertidos novamente para segundos.

Por exemplo:
```
delay_index = [1 3 5 8]
Ts = 1 ms
```

resulta em:
```
delay = [0 2 4 7] ms
```

Então: $\Large \space = \space (n_l \space - \space 1)T_s$

### 3.2 Segunda etapa: geração dos ganhos
Agora o programa calcula os ganhos dos caminhos.

Primeiro:
```
alpha = log(10^(gain.attenuation/10))/setup.delayspread;
```

Esse parâmetro `alpha` controla a velocidade com que a potência média diminui com o atraso.

A ideia é: $\Large P_l \space \varpropto \space e^{-\alpha \tau_l}$
​
Quanto maior o atraso, menor tende a ser a potência média do caminho.

#### Variância do ganho
Depois:
```
gain_variance = exp(-alpha*delay);
```

Isso gera uma variância diferente para cada caminho.
Por exemplo:
```
delay             gain_variance
--------------------------------
0 ms                 1.00
2 ms                 0.75
4 ms                 0.50
6 ms                 0.30
```

Os valores acima são apenas ilustrativos.

A ideia é: $\Large \sigma_l^2 \space = \space e^{-\alpha \tau_l}$
Portanto:
> **Quanto maior o atraso, menor tende a ser a potência do caminho.**

#### Distribuição Rayleigh para os ganhos
Depois:
```
gain_tap = raylrnd(sqrt(gain_variance*2/(4 - pi)));
```

Aqui o código gera os ganhos utilizando uma **distribuição Rayleigh**.

Isso é bastante comum em modelos de canais sem uma componente dominante de linha de visada.

Assim, cada caminho recebe um ganho aleatório:$$\Large h_l​ \space \thicksim \space Rayleigh(\sigma_l)$$
Por exemplo:
```
gain_tap =
    0.52
    0.31
    0.18
    0.11
    0.04
```

Novamente, valores ilustrativos.

#### Normalização dos ganhos
Depois:
```
gain_tap = gain_tap/norm(gain_tap);
```
`norm(gain_tap)` calcula:$$\Large \sqrt{∣h_1​∣^2+∣h_2​∣^2+ \dots +∣h_L​∣^2​}$$
Então o vetor é normalizado:$$\Large h \space \leftarrow \space \frac{h}{||h||}$$

Isso faz com que: $\large \sum_l ​|h_l​|^2\space = \space 1$
aproximadamente.

Portanto, os ganhos são escalados para que a energia total dos taps seja unitária.

### 3.3 Terceira etapa: Doppler
Agora chegamos ao Doppler.

O código define:
```
c = 1500;
```

Esse é o valor aproximado da velocidade do som na água:$$\large 
c \approx 1500 \space m/s$$
O fator Doppler é relacionado à velocidade relativa: $\large a \space = \space \frac{v}{c}$
onde:
- $\Large v$ = velocidade relativa;
- $\Large c$ = velocidade do som.

#### Caso `doppler.type = 'none'`
```
if(strcmp(doppler.type,'none'))
    a_max = 0;
    varargout{1} = 1;
    varargout{2} = 1;
```

Nesse caso:
$\large a_{max}​ = 0$
Portanto:
$\large 1+a_{max}​=1$
e não há alteração nos atrasos.

Também temos:
```
Q = 1
M = 1
```

#### Caso Doppler uniforme
```
elseif(strcmp(doppler.type,'uniform'))
```

Primeiro:
```
a_max = doppler.velocity/c;
```

Se:
```
doppler.velocity = 15;
```

então:
$\large a_{max​}=\frac{15}{1500}​=0.01$
Ou seja:
$\large 1+a_{max}​=1.01$
Depois:
```
[Q,M] = rat(1+a_max);
```
`rat` procura uma aproximação racional.

Por exemplo:
$\large 1.01 \space = \space \frac{101}{100}​$
Então:
```
Q = 101
M = 100
```
Esses valores são utilizados para representar o fator de amostragem associado ao Doppler.

#### Caso Doppler não uniforme
A intenção do código é:
```
elseif (strcmp(type,'non-uniform'))
```

calcular:
```
a = doppler.velocity/c;
a_max = max(a);
```

Porém existe um **erro importante**:
```
strcmp(type,'non-uniform')
```

deveria provavelmente ser:
```
strcmp(doppler.type,'non-uniform')
```

Porque `type` não foi definido anteriormente.

Então, do jeito que está, essa condição provavelmente causará erro.

#### Correção dessa parte
Deveria ser:
```
elseif(strcmp(doppler.type,'non-uniform'))
    a = doppler.velocity/c;
    a_max = max(a);
    [Q,M] = rat(1+a_max);
    varargout{1} = Q;
    varargout{2} = M;
```

Aqui cada caminho pode possuir uma velocidade diferente.

Por exemplo:
```
doppler.velocity = [5 10 15 20];
```

Então:
$\large a \space = \space \Big{[} \frac{5}{1500}, \frac{10}{1500}, \frac{15}{1500}, \frac{20}{1500} \Big{]}$
e o maior fator é usado para definir `a_max`.

**Entretanto**, há outra questão: apesar do nome `non-uniform`, o código usa `a_max` para corrigir todos os atrasos. Portanto, ele não está aplicando um fator Doppler individual a cada caminho nessa parte.
#### Correção dos atrasos pelo Doppler
Depois vem:
```
delay_index_bar = ceil(delay_index./(1 + a_max));
```

Esse é o ponto em que o Doppler modifica os atrasos.
Sem Doppler:
$$\large a_{max}​=0$$
então:
$\large \bar{n}_l \space = \space \Big{[} \frac{n_l}{1} \Big{]} = n_l$​

Com Doppler:
$\large \bar{n}_l \space = \space \Big{[} \frac{n_l}{1+a_{max}} \Big{]}$

#### Convertendo os novos índices para tempo
Depois:
```
delay_bar = (delay_index_bar - 1)*setup.Ts;
```

Então:
$\Large \bar{\tau}_l \space = \space (\bar{n}_l - 1) T_s$
​

`delay_bar` é, portanto, o vetor de atrasos **depois da correção relacionada ao Doppler**.

#### Construção de `h`
Finalmente:
```
h(delay_index_bar) = gain_tap;
```

Essa linha é extremamente importante.
Ela cria a resposta ao impulso do canal.

Imagine:
```
delay_index_bar = [1 4 7]
gain_tap        = [0.8 0.5 0.2]
```

Então:
```
h =
[0.8  0  0  0.5  0  0  0.2]
```

Graficamente:
```
amplitude
   |
0.8| *
   |
0.5|       *
   |
0.2|             *
   |
   +------------------------> tempo
     0   3   6 ...
```

Isso representa os três caminhos do canal.

## 4. Fluxo completo
Podemos resumir o programa assim:
```
                 ENTRADAS
                    |
       +------------+------------+
       |            |            |
     setup         gain        delay
       |            |            |
       +------------+------------+
                    |
                    v
        Geração dos intervalos
             entre caminhos
                    |
                    v
       Distribuição exponencial
                    |
                    v
          Quantização temporal
                    |
                    v
          Atrasos acumulados
                    |
                    v
          delay_index / delay
                    |
                    v
       Cálculo da atenuação média
                    |
                    v
        Variância dos caminhos
                    |
                    v
         Ganhos Rayleigh
                    |
                    v
         Normalização dos ganhos
                    |
                    v
             DOPPLER
                    |
          +---------+---------+
          |         |         |
        none     uniform   non-uniform
          |         |         |
          +---------+---------+
                    |
                    v
          Correção dos atrasos
                    |
                    v
               delay_bar
                    |
                    v
       Colocação dos ganhos em h
                    |
                    v
                  SAÍDA
```

## 5. O que cada saída representa?
## `h`
É a **resposta ao impulso discreta do canal**.
Exemplo:
```MATLAB
h = 0.82    0    0.51    0    0    0.23
```

Isso significa que existem três caminhos significativos:
```
amostra 1 → ganho 0.82
amostra 3 → ganho 0.51
amostra 6 → ganho 0.23
```

## `delay_bar`

É o atraso de cada caminho depois da consideração do fator Doppler.
Exemplo
```
delay_bar =
    0
    0.002
    0.005
```

Ou seja:
$ \large \begin{aligned} \tau_1​=0 \\ \tau_2​=2ms \\ \tau_3​=5ms \end{aligned}$
## `gain_tap`

É o ganho de cada caminho:
```
gain_tap =
    0.82
    0.51
    0.23
```
Existe uma correspondência:
```
delay_bar(1) <-> gain_tap(1)
delay_bar(2) <-> gain_tap(2)
delay_bar(3) <-> gain_tap(3)
```

## `varargout{1}`

É:
```
Q
```
o fator de down sampling associado ao fator Doppler racionalizado.

## `varargout{2}`

É:
```
M
```
o fator de upsampling.

Eles aparecem porque o fator:
$\Large 1+a_{max}$​
é aproximado por:
$\Large \frac{Q}{M}$

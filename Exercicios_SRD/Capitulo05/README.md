## Modulação por amplitude
Considerando uma onda portadora $\large c(t)$ definida por
$$
\large c(t) = A_c \space cos(2 \pi f_c t)
$$

onde $\large A_c$ é a amplitude da portadora e $\large f_c$ e a frequência da portadora. $\large m(t)$ é a mensagem a ser transmitida.
O sinal modulado (Amplitude Modulated-AM) pode ser definido por meia da seguinte função
$$
\large s(t) = A_c[1+k_a m(t)] \space cos(2\pi f_c t)
$$
- Ac​ **(Amplitude da Portadora):** Representa a **amplitude máxima da onda portadora** senoidal. Na equação que define o sinal AM, Ac​ atua como um fator de escala para o sinal transmitido.
- ka​ **(Sensibilidade de Amplitude):** É a constante que define a **sensibilidade do modulador** encarregado de gerar o sinal AM. Ela determina o quanto a amplitude da portadora será alterada em resposta às variações do sinal de mensagem m(t). A unidade de ka​ é o inverso da unidade do sinal de mensagem (geralmente $\large {volts}^{-1}$ se a mensagem for em Volts)

Normalmente a Amplitude da portadora e o sinal da mensagem são medidos em $\large volts$ e $\large k_a$ é medido em $\large {volts}^{-1}$.
**Condições de Operação**
Para que o sinal seja transmitido e recuperado corretamente sem distorção (usando, por exemplo, um detector de envelope), as fontes destacam que:
1. **Evitar Sobremodulação:** O valor absoluto de ka​m(t) deve ser mantido **menor que a unidade** para todo $\large t \space (∣ka​m(t)∣<1)$. Se essa condição for violada, o envelope do sinal pode passar por zero e tornar-se negativo, resultando em uma portadora "reversora de fase" e distorção no sinal recuperado.
2. A frequência da portadora $\large f_c$ deve ser muito maior que a frequência da mensagem $\large W$, ou seja $$\large f_c >> W$$

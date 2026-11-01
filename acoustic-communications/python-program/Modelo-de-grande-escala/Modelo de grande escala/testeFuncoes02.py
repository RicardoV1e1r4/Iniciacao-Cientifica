import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
f1 = 1000  # Frequência em kHz
d1 = 1      # Distância em km
Ks = 1      # Coeficiente de absorção
delta_d = 0.1  # Incremento de distância em km
x_start = 100; x_end = 10000
y_start = 0; y_end = 0
#------------------------------------------------

# Parâmetros variáveis
#------------------------------------------------
f_tl = np.array([ 0.1, 1, 5, 10, 100 ])
d = np.arange(100, 100000, 0.1)  # Distâncias de 0.1 km a 10 km
tempo = np.arange(0, 20, 0.1)  # Tempo de 0 a 20 segundos
#------------------------------------------------


#RAZÃO DA PERDA DE PROPAGAÇÃO
#------------------------------------------------
xposition = 0.5 * np.cos(tempo)   # Posição x variando com o tempo
distance0, distance1, variacao = uwf.distanceVariation(x_start, x_end, xposition, y_start, y_end, 0)
# print(distance0, distance1)

beta = (1 + (variacao / distance0)) ** Ks * (10 ** ((uwf.ThorpsFormula(f1) / 10000) * variacao))

fig, axes = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=False, sharex=False, sharey=False)

# Agora você acessa como uma matriz [linha, coluna]
ax1 = axes[0]
ax2 = axes[1]

ax1.plot(tempo, beta)
ax1.set_xlabel('Tempo (s)')
ax1.set_ylabel('Beta')
ax1.set_title('Variação da Perda de Transmissão')
ax1.grid(True)

ax2.plot(tempo, 1 / (beta**0.5))
ax2.set_xlabel('Tempo (s)')
ax2.set_ylabel('A_LOS')
ax2.set_title('Variação da Perda de Transmissão')
ax2.grid(True)

plt.tight_layout()
plt.show()

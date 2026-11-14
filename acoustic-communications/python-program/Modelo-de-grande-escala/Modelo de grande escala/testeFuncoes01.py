import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
f1 = 10     # Frequência em kHz
Ks = 2      # Coeficiente de absorção
L0 = 1      # Perda de referência (escala linear)
#------------------------------------------------

# Vetores
#------------------------------------------------
f       = np.array([ 0.1, 1, 2, 5, 10, 20, 50, 100 ])   # Frequências de 0.1 kHz a 100 kHz
f_tl    = np.array([ 0.1, 1, 5, 10, 20, 50, 100 ])      # Frequências de 0.1 kHz a 100 kHz
d       = np.arange(0.1, 100, 0.1)                      # Distâncias de 0.1 km a 100 km

#------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=False, sharex=False, sharey=False)

# Agora você acessa como uma matriz [linha, coluna]
ax1 = axes[0]
ax2 = axes[1]

# PERDA DE PROPAGAÇÃO
#------------------------------------------------
ax1.set_xscale('log')
ax1.set_yscale('log')

for i in f_tl:
    loss = uwf.conventionalPropagationLoss(i, d)
    ax1.plot(d * 1000, loss, label=f'{i} kHz')

ax1.set_xlabel('Distância (m)')
ax1.set_ylabel('Perda de Propagação (dB)')
ax1.set_title('Perda de Propagação Convencional')
ax1.set_xlim(100, 100000)
ax1.set_ylim(40, 120)
ax1.legend()
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

# PERDA DE TRANSMISSÃO
#------------------------------------------------
ax2.set_xscale('log')
ax2.set_yscale('log')

for i in f_tl:
    loss = 10 * np.log10(uwf.transmissionLoss(i, d, L0, Ks))
    ax2.plot(d * 1000, loss, label=f'{i} kHz')

ax2.set_xlabel('Distância (m)')
ax2.set_ylabel('Perda de Transmissão (dB)')
ax2.set_title('Perda de Transmissão')
ax2.set_xlim(100, 100000)
ax2.set_ylim(30, 150)
ax2.legend()
ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
f1      = 10    # Frequência em kHz
L0      = 1     # Perda de referência (escala linear)
Ks      = 1.5   # Fator de espalhamento (1 = cilíndrico, 2 = esférico, 1.5 = prático)
v_R     = 0     # Velocidade do receptor em m/s
v_T     = 0     # Velocidade do transmissor em m/s
theta_R = 0     # Direção do receptor em graus (0 = +x, 90 = +y, 180 = -x)
theta_T = 0     # Direção do transmissor em graus
#------------------------------------------------

# Componentes da velocidade do receptor, em m/s
v_Rx = v_R * np.cos(np.deg2rad(theta_R))
v_Ry = v_R * np.sin(np.deg2rad(theta_R))

# Componentes da velocidade do transmissor, em m/s
v_Tx = v_T * np.cos(np.deg2rad(theta_T))
v_Ty = v_T * np.sin(np.deg2rad(theta_T))

# Posição inicial do receptor, em km (as seguintes saem de v_R, theta_R e tempo)
x_R0 = 0
y_R0 = 0

# Posição inicial do transmissor, em km (as seguintes saem de v_T, theta_T e tempo)
x_T0 = 0.5
y_T0 = 0.5

# Vetores
#------------------------------------------------
f_tl    = np.array([ 0.1, 1, 5, 10, 100 ])  # Frequências de 0.1 kHz a 100 kHz
d       = np.arange(0.1, 100, 0.1)          # Distâncias de 0.1 km a 10 km
tempo   = np.arange(0, 60, 0.1)             # Tempo de 0 a 20 segundos
#------------------------------------------------


#RAZÃO DA PERDA DE PROPAGAÇÃO
#------------------------------------------------
def distancia(t):
    """Distância transmissor-receptor em km no instante t (em s)."""
    x_R1 = x_R0 + v_Rx * t / 1000
    y_R1 = y_R0 + v_Ry * t / 1000
    x_T1 = x_T0 + v_Tx * t / 1000
    y_T1 = y_T0 + v_Ty * t / 1000
    return np.hypot(x_T1 - x_R1, y_T1 - y_R1)


TL = (uwf.transmissionLoss(f1, distancia(tempo), 1, Ks))

fig, axes = plt.subplots(1, 2, figsize=(10, 6), constrained_layout=False, sharex=False, sharey=False)

# Agora você acessa como uma matriz [linha, coluna]
ax1 = axes[0]
ax2 = axes[1]

# ax1.set_yscale('log')
# ax1.set_xscale('log')

# ax2.set_yscale('log')
# ax2.set_xscale('log')

ax1.plot(d * 1000, TL)
ax1.set_xlabel('Distância (m)')
ax1.set_ylabel('TL(d)')
ax1.set_title('Variação da Perda de Transmissão')
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

ax2.plot(d * 1000, 1 / TL)
ax2.set_xlabel('Distância (m)')
ax2.set_ylabel('1/TL(d)')
ax2.set_title('Variação da Perda de Transmissão')
ax2.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.show()

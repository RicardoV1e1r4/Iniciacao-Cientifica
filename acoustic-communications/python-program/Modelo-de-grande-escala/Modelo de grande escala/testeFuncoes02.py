import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
f1      = 10    # Frequência em kHz
Ks      = 1.5   # Fator de espalhamento (1 = cilíndrico, 2 = esférico, 1.5 = prático)
v_R     = 5     # Velocidade do receptor em m/s
v_T     = 0     # Velocidade do transmissor em m/s
theta_R = 0     # Direção do receptor em graus (0 = +x, 90 = +y, 180 = -x)
theta_T = 0     # Direção do transmissor em graus
Tbl     = 60    # Duração do bloco em s
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
x_T0 = 0.4
y_T0 = 0.4

# Vetores
#------------------------------------------------
f_tl = np.array([ 0.1, 1, 5, 10, 100 ])     # Frequências de 0.1 kHz a 100 kHz
d = np.arange(0.1, 100, 0.1)                # Distâncias de 0.1 km a 100 km
tempo = np.arange(0, Tbl, 0.1)               # Tempo de 0 a Tbl segundos
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


distance0 = distancia(0)            # Distância no início
distance1 = distancia(tempo)        # Distância ao longo do tempo
dr = distance1 - distance0          # Variação da distância em relação ao início

# Razão TL(d)/TL(d0): alpha em dB/km e dr em km, por isso o expoente é alpha*dr/10
beta = (1 + (dr / distance0)) ** Ks * (10 ** ((uwf.ThorpsFormula(f1) / 10) * dr))

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

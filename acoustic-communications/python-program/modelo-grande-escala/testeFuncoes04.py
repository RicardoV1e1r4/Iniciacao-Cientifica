import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
f1      = 10    # Frequência em kHz
L0      = 1     # Perda de referência (escala linear)
Ks      = 2     # Fator de espalhamento (1 = cilíndrico, 2 = esférico, 1.5 = prático)
v_R     = 10     # Velocidade do receptor em m/s
v_T     = 0     # Velocidade do transmissor em m/s
theta_R = 0     # Direção do receptor em graus (0 = +x, 90 = +y, 180 = -x)
theta_T = 0     # Direção do transmissor em graus
T_max   = 10    # Tempo máximo da simulação em segundos
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
x_T0 = 6.7
y_T0 = 1.8

# Vetores
#------------------------------------------------
d       = np.logspace(np.log10(0.1), np.log10(100), 120)   # Distâncias de 0.1 km a 100 km (malha da superfície)
tempo   = np.arange(0, T_max, 0.1)                            # Tempo de 0 a 60 segundos (trajetória)
t_surf  = np.linspace(tempo[0], tempo[-1], 60)             # Tempo mais esparso para a superfície
#------------------------------------------------

D, T = np.meshgrid(d, t_surf, indexing='ij')  # Malha (linhas = distância, colunas = tempo)


# PERDA DE TRANSMISSÃO (dB)
#------------------------------------------------
def transmissionLoss_dB(f, d_km, L0, Ks):
    """TL em dB: 10*log10(L0) + Ks*10*log10(d) + alpha*d.

    Equivale a 10*log10(uwf.transmissionLoss(f, d, L0, Ks)), mas evita overflow
    do valor linear em distâncias grandes.
    """
    return 10 * np.log10(L0) + Ks * 10 * np.log10(d_km * 1000) + uwf.ThorpsFormula(f) * d_km


def distancia(t):
    """Distância transmissor-receptor em km no instante t (em s)."""
    x_R1 = x_R0 + v_Rx * t / 1000
    y_R1 = y_R0 + v_Ry * t / 1000
    x_T1 = x_T0 + v_Tx * t / 1000
    y_T1 = y_T0 + v_Ty * t / 1000
    return np.hypot(x_T1 - x_R1, y_T1 - y_R1)


# Superfície: TL em função da distância, repetida ao longo do tempo
TL_surf = transmissionLoss_dB(f1, D, L0, Ks)

upper_limit = 130
lower_limit = 40
TL_surf_recortado = np.ma.masked_outside(TL_surf, lower_limit, upper_limit)

# Trajetória: distância real entre transmissor e receptor ao longo do tempo
d_t  = distancia(tempo)
TL_t = transmissionLoss_dB(f1, d_t, L0, Ks)

# Propagação da frente de onda (pulso emitido em t = 0)
# Após t_p segundos, a frente de onda está em r = c*t_p e acumulou TL(r).
# Começa em r = 1 m (distância de referência, onde TL = 10*log10(L0)).
c    = 1500                                   # Velocidade do som na água em m/s
d0   = distancia(0)                           # Distância do enlace no instante da emissão, em km
tau  = d0 * 1000 / c                          # Tempo de propagação até o receptor, em s
t_p  = np.linspace(1 / c, tau, 500)           # Tempo de propagação de 1 m até o receptor
r_p  = c * t_p / 1000                         # Posição da frente de onda em km
TL_p = transmissionLoss_dB(f1, r_p, L0, Ks)


# ESCALA LOGARÍTMICA
#------------------------------------------------
# Os eixos 3D do matplotlib não suportam set_xscale('log') de forma confiável.
# Por isso a distância é plotada em log10 e os rótulos mostram o valor real.
def log_axis(axis, vmin, vmax, minor=True):
    """Configura um eixo (já em log10) com marcações nas potências de 10."""
    decades = np.arange(np.floor(np.log10(vmin)), np.ceil(np.log10(vmax)) + 1)
    ticks = decades
    if minor:
        # Adiciona marcações intermediárias (2 e 5 de cada década)
        extra = np.log10(np.concatenate([[2 * 10**k, 5 * 10**k] for k in decades]))
        ticks = np.sort(np.concatenate([decades, extra]))
    ticks = ticks[(ticks >= np.log10(vmin) - 1e-9) & (ticks <= np.log10(vmax) + 1e-9)]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f'{10**v:g}'))


logD  = np.log10(D * 1000)     # distância em metros
logdt = np.log10(d_t * 1000)   # distância da trajetória em metros

# GRÁFICOS
#------------------------------------------------
fig = plt.figure(figsize=(20, 7))
fig.suptitle(f'Perda × Distância × Tempo (f = {f1} kHz, Ks = {Ks})', fontsize=14)

ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1, 2, 2)
# ax3 = fig.add_subplot(1, 3, 3)

# --- 3D: superfície TL(d, t) + trajetória do enlace ---
surf1 = ax1.plot_surface(logD, T, TL_surf_recortado, cmap='plasma', alpha=0.5,
                         linewidth=0, antialiased=True)
ax1.plot(logdt, tempo, TL_t, color='red', linewidth=2.5, label='Trajetória do enlace')

log_axis(ax1.xaxis, d[0] * 1000, d[-1] * 1000, minor=False)
ax1.set_xlabel('Distância (m)')
ax1.set_ylabel('Tempo (s)')
ax1.set_zlabel('TL (dB)')
ax1.set_title('Perda de Transmissão')
ax1.legend(loc='upper left', fontsize=8)
fig.colorbar(surf1, ax=ax1, shrink=0.6, pad=0.1, label='TL (dB)')

# --- 2D: TL ao longo do tempo para a trajetória ---
ax2.plot(tempo, TL_t, color='red')
ax2.set_xlabel('Tempo (s)')
# ax2.set_ylabel('TL (dB)')
ax2.set_title('Perda de Transmissão ao longo do tempo')
ax2.grid(True)

# --- 2D: perda acumulada pela frente de onda até chegar ao receptor ---
fig2, ax3 = plt.subplots(figsize=(8, 5))
ax3.plot(t_p, TL_p, color='blue', label='Frente de onda')
ax3.plot(tau, TL_p[-1], 'o', color='red',
         label=f'Receptor: τ = {tau:.3f} s, TL = {TL_p[-1]:.1f} dB')
ax3.axvline(tau, color='red', linestyle='--', linewidth=1)
ax3.set_xlabel('Tempo de propagação (s)')
ax3.set_ylabel('TL acumulada (dB)')
ax3.set_title(f'Propagação do pulso emitido em t = 0 (c = {c} m/s)')
ax3.legend(loc='lower right', fontsize=8)
ax3.grid(True)

# Eixo superior com a posição da frente de onda em metros
ax3_top = ax3.secondary_xaxis('top', functions=(lambda t: t * c, lambda r: r / c))
ax3_top.set_xlabel('Posição da frente de onda (m)')

plt.tight_layout()
plt.show()

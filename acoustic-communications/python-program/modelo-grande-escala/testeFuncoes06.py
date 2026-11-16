import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
Ks = 2      # Fator de espalhamento (1 = cilíndrico, 1.5 = prático, 2 = esférico)
L0 = 1      # Perda de referência (escala linear)
#------------------------------------------------

# Parâmetros variáveis (espaçamento logarítmico, como nos eixos do testeFuncoes01)
#------------------------------------------------
f       = np.logspace(np.log10(0.1), np.log10(100), 120)    # Frequências de 0.1 kHz a 100 kHz
f_tl    = np.logspace(np.log10(0.1), np.log10(60), 8)      # Frequências de teste para perda de propagação
d       = np.logspace(np.log10(0.1), np.log10(100), 120)    # Distâncias de 0.1 km a 100 km
#------------------------------------------------

# Malha (linhas = frequência, colunas = distância)
F, D = np.meshgrid(f, d, indexing='ij')

# PERDA DE PROPAGAÇÃO (dB)
#------------------------------------------------
PL = uwf.conventionalPropagationLoss(F, D)

PL_recortado = np.where((PL >= 40) & (PL <= 130), PL, np.nan)

# PERDA DE TRANSMISSÃO (dB)
#------------------------------------------------
# Calculada em dB dentro da função (evita overflow do valor linear).
TL = uwf.transmissionLoss(F, D, L0, Ks)

# ESCALA LOGARÍTMICA
#------------------------------------------------
# Os eixos 3D do matplotlib não suportam set_xscale('log') de forma confiável.
# Por isso os dados são plotados em log10 e os rótulos dos eixos mostram o
# valor real (equivalente ao type='log' do plotly no testeFuncoes05).
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


logD  = np.log10(D * 1000)   # distância em metros
logF  = np.log10(F)          # frequência em kHz
logd  = np.log10(d * 1000)

# GRÁFICOS 3D
#------------------------------------------------
fig = plt.figure(figsize=(16, 7.5))
fig.suptitle('Perda × Distância × Frequência', fontsize=14)

ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1, 2, 2, projection='3d')

# --- Perda de propagação convencional (recortada entre 40 e 130 dB) ---
# plot_surface não lida bem com NaN; usamos um array mascarado para ocultar
# as faces fora do intervalo.
PL_mask = np.ma.masked_invalid(np.log10(PL_recortado))
surf1 = ax1.plot_surface(logD, logF, PL_mask, cmap='plasma', alpha=0.5,
                         vmin=np.log10(40), vmax=np.log10(130),
                         linewidth=0, antialiased=True)

# Curvas para frequências fixas, projetadas na parede de frequência máxima
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'magenta', 'gray']
y_parede = np.log10(f[-1])
for fi, cor in zip(f_tl, colors):
    PL1 = uwf.conventionalPropagationLoss(fi, d)
    PL1_recortado = np.where((PL1 >= 40) & (PL1 <= 130), PL1, np.nan)
    ax1.plot(logd, np.full_like(d, y_parede), np.log10(PL1_recortado),
             color=cor, linewidth=2.5, label=f'f = {fi:.2f} kHz')

ax1.set_title('Perda de Propagação Convencional')
ax1.legend(loc='upper left', fontsize=8)

# --- Perda de transmissão ---
surf2 = ax2.plot_surface(logD, logF, np.log10(TL), cmap='plasma',
                         linewidth=0, antialiased=True)
ax2.set_title(f'Perda de Transmissão (Ks = {Ks})')

# --- Eixos, escalas e câmera (iguais nos dois gráficos) ---
for ax, z in ((ax1, PL_recortado), (ax2, TL)):
    ax.set_xlabel('Distância (m)')
    ax.set_ylabel('Frequência (kHz)')
    ax.set_zlabel('Perda (dB)')

    ax.set_xlim(np.log10(d[0] * 1000), np.log10(d[-1] * 1000))
    ax.set_ylim(np.log10(f[0]), np.log10(f[-1]))
    zmin, zmax = np.nanmin(z), np.nanmax(z)
    ax.set_zlim(np.log10(zmin), np.log10(zmax))

    log_axis(ax.xaxis, d[0] * 1000, d[-1] * 1000, minor=False)
    log_axis(ax.yaxis, f[0], f[-1], minor=False)
    log_axis(ax.zaxis, zmin, zmax)

    # Câmera equivalente a eye=(-1.6, -1.6, 0.9) do plotly
    ax.view_init(elev=22, azim=-135)
    ax.set_box_aspect(None, zoom=0.88)   # evita que os rótulos dos eixos sejam cortados

# Barras de cores com os valores reais em dB
for surf, ax, label in ((surf1, ax1, 'PL (dB)'), (surf2, ax2, 'TL (dB)')):
    cb = fig.colorbar(surf, ax=ax, shrink=0.7, pad=0.1)
    cb.set_label(label)
    vmin, vmax = surf.get_clim()
    log_axis(cb.ax.yaxis, 10**vmin, 10**vmax)
    cb.solids.set_alpha(1)

plt.tight_layout()
plt.show()

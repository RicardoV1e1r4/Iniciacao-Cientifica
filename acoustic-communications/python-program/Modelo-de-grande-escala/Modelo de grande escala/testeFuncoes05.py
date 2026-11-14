import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import underwaterFunctions as uwf


# Parâmetros estaticos
#-------------------------------------------------
Ks = 2      # Fator de espalhamento (1 = cilíndrico, 1.5 = prático, 2 = esférico)
L0 = 1      # Perda de referência (escala linear)
#------------------------------------------------

# Parâmetros variáveis (espaçamento logarítmico, como nos eixos do testeFuncoes01)
#------------------------------------------------
f       = np.logspace(np.log10(0.1), np.log10(100), 120)    # Frequências de 0.1 kHz a 100 kHz
f_tl    = np.logspace(np.log10(0.1), np.log10(70), 8)       # Frequências de teste para perda de transmissão
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
# Calculada diretamente em dB: 10*log10(L0) + Ks*10*log10(d) + alpha*d.
# Equivale a 10*log10(uwf.transmissionLoss(...)), mas evita overflow do valor
# linear em frequências e distâncias altas (ex.: 100 kHz a 100 km -> ~10^350).
# TL = L0 * ((d * 1000) ** Ks) * (10**((alpha / 10) * d))

TL = 10 * np.log10(L0) + Ks * 10 * np.log10(D * 1000) + uwf.ThorpsFormula(F) * D

# TL1 = uwf.conventionalPropagationLoss(f_tl, d)

# GRÁFICOS 3D
#------------------------------------------------
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'scene'}]],
    subplot_titles=('Perda de Propagação Convencional', f'Perda de Transmissão (Ks = {Ks})'),
    horizontal_spacing=0.02,
)

hover = ('Distância:    %{x:.1f} m<br>'
         'Frequência:   %{y:.2f} kHz<br>'
         'Perda:        %{z:.1f} dB<extra></extra>')

fig.add_trace(
    go.Surface(x=D * 1000, y=F, z=PL_recortado, colorscale='plasma', opacity=0.5,
               colorbar=dict(title='PL (dB)', x=0.45, len=0.8),
               hovertemplate=hover),
    row=1, col=1,
)

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'magenta', 'gray']
for i in f_tl:
    PL1 = uwf.conventionalPropagationLoss(i, d)
    PL1_recortado = np.where((PL1 >= 40) & (PL1 <= 130), PL1, np.nan)

    fig.add_trace(
        go.Scatter3d(x=d * 1000, y=np.full_like(d, d[-1]), z=PL1_recortado, mode='lines',
                    line=dict(color=colors[f_tl.tolist().index(i)], width=5),
                    name=f'f = {i} kHz'),
        row=1, col=1,
    )

fig.add_trace(
    go.Surface(x=D * 1000, y=F, z=TL, colorscale='plasma',
               colorbar=dict(title='TL (dB)', x=1.0, len=0.8),
               hovertemplate=hover),
    row=1, col=2,
)

scene_layout = dict(
    xaxis=dict(title='Distância (m)', type='log'),
    yaxis=dict(title='Frequência (kHz)', type='log'),
    zaxis=dict(title='Perda (dB)', type='log'),
    camera=dict(eye=dict(x=-1.6, y=-1.6, z=0.9)),
)

fig.update_layout(
    title='Perda × Distância × Frequência',
    scene=scene_layout,
    scene2=scene_layout,
    width=1400, height=700,
)

fig.show()

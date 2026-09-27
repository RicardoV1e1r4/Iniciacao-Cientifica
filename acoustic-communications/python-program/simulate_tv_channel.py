import os
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

#os.makedirs('/workspace/scratch', exist_ok=True)
sns.set_theme(style='white', palette='colorblind', font='DejaVu Sans')

# Parameters
# Absolute time t from 0 to 1 second
t = np.linspace(0, 1.0, 200) # [s]
# Delay tau from 0 to 30 ms
tau_ms = np.linspace(0, 30, 300) # [ms]

T, TAU = np.meshgrid(t, tau_ms)

# Multipath parameters according to Eq 1.14:
# h(t; \tau) = \sum A_p \delta(\tau - (\tau_p - a_p t))
# We model delta as a Gaussian pulse with pulse width sigma_tau = 0.4 ms
sigma = 0.4 # [ms]

# Paths configuration (A_p, tau_p in ms, a_p Doppler scaling factor)
paths = [
    {'A': 1.0,  'tau0': 5.0,  'a': 0.0008, 'name': 'Raio Direto'},
    {'A': -0.75, 'tau0': 11.0, 'a': -0.0015, 'name': 'Reflexão Superfície'},
    {'A': 0.5,  'tau0': 18.0, 'a': 0.0005, 'name': 'Reflexão Fundo'},
    {'A': -0.35, 'tau0': 24.0, 'a': -0.0022, 'name': 'Reflexão Dupla (Superfície-Fundo)'}
]

H = np.zeros_like(T)

for p in paths:
    # tau_p(t) = tau0 - a_p * t * 1000 (convert t to ms scaling)
    # Note: tau = tau_p - a_p * t => tau_p(t) in ms: tau0 - a_p * t * 1000
    tau_p_t = p['tau0'] - p['a'] * T * 1000.0
    pulse = p['A'] * np.exp(-0.5 * ((TAU - tau_p_t) / sigma)**2)
    H += pulse

# Create 3D visualization
fig = plt.figure(figsize=(14, 7), dpi=150)

# Subplot 1: 3D Surface
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax1.plot_surface(T, TAU, np.abs(H), cmap='viridis', edgecolor='none', alpha=0.9)
ax1.set_title('Resposta ao Impulso $h(t, \\tau)$ (Superfície 3D)\nAtrasos Variantes com o Tempo Absoluto', fontsize=12, fontweight='bold', pad=12)
ax1.set_xlabel('Tempo Absoluto $t$ [s]', fontsize=10, labelpad=8)
ax1.set_ylabel('Atraso $\\tau$ [ms]', fontsize=10, labelpad=8)
ax1.set_zlabel('Magnitude $|h(t, \\tau)|$', fontsize=10, labelpad=8)
ax1.view_init(elev=28, azim=-55)

# Subplot 2: 2D Intensity Map / Waterfall Projection
ax2 = fig.add_subplot(1, 2, 2)
contour = ax2.pcolormesh(T, TAU, np.abs(H), cmap='viridis', shading='auto')
cbar = fig.colorbar(contour, ax=ax2, pad=0.03)
cbar.set_label('Magnitude $|h(t, \\tau)|$', fontsize=10)

ax2.set_title('Visão Superior (Mapa de Intensidade Tempo-Atraso)\nTrajetórias Doppler dos Percursos Múltiplos', fontsize=12, fontweight='bold', pad=12)
ax2.set_xlabel('Tempo Absoluto $t$ [s]', fontsize=10)
ax2.set_ylabel('Atraso $\\tau$ [ms]', fontsize=10)

# Annotate trajectories on 2D map
for p in paths:
    tau_line = p['tau0'] - p['a'] * t * 1000.0
    ax2.plot(t, tau_line, '--', color='white', linewidth=1.2, alpha=0.7)

plt.tight_layout(pad=2.0)

#plt.show()

# Save to scratch
scratch_path = 'D:/Documentos/GitHub/Cognitive-Radios/acoustic-communications/python-program/time_varying_impulse_response_3d.png'
fig.savefig(scratch_path, bbox_inches='tight', dpi=150)
plt.close(fig)

print(f"Generated successfully at {scratch_path}")

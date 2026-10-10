# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:15:08 2026

@author: Ricardo Alexandre
"""

import matplotlib.pyplot as plt
import numpy as np
import acousticChannel2 as ac

fs = 1000  # Hz (Taxa de amostragem)
Tbl = 0.025

setup = {
    "Ts": 1 / fs,
    "paths": 10,
    "delayspread": 0.025
}

gain = {
    "mode": "deterministic",  # Ativa o modelo determinístico por Taylor
    # "mode": "stochastics",  # Ativa o modelo determinístico por Taylor
    "attenuation": 20,
    "dt": 50 / 48000,          # Intervalo dt temporal (~1.04 ms)
    "T_bl": 0.025              # Duração do bloco OFDM (25 ms)
}

delay = {
    "mean": 0.003
}

doppler = {
    "type": "uniform",
    "velocity": 1.5,
    "f_d": 2.0                 # Espalhamento Doppler
}

h, delay_bar, gain_tap, Q, M = ac.acoustic_channel(setup, gain, delay, doppler)

print("\nh =", h)
print("\ndelay_bar =", delay_bar)
print("\ngain_tap =", gain_tap)
print("\nQ =", Q)
print("\nM =", M)
print("\nEnergia total dos taps:", np.dot(gain_tap, gain_tap))

# Plotagem dos resultados
fig, axes = plt.subplots(1, 2, figsize=(10, 6))

axes[0].stem(h)
axes[0].set_xlabel("Índice da amostra n")
axes[0].set_ylabel("Amplitude")
axes[0].set_title("Resposta impulsiva discreta do canal")
axes[0].grid(True)

axes[1].stem(delay_bar, gain_tap)
axes[1].set_xlabel("Atraso (s)")
axes[1].set_ylabel("Amplitude (ganho)")
axes[1].set_title("Resposta do canal multipercurso (Taylor)")
axes[1].grid(True)

plt.tight_layout()
plt.show()
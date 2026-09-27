# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:15:08 2026

@author: Ricardo Alexandre
"""

import matplotlib.pyplot as plt
import acousticChannel as ac

setup = {
    "Ts": 0.001,
    "paths": 8,
    "delayspread": 0.01}

gain = {
    "attenuation": 15}

delay = {
    "mean": 0.003}

doppler = {
    "type": "uniform",
    "velocity": 15}

h, delay_bar, gain_tap, Q, M = ac.acoustic_channel(setup, gain, delay, doppler)

print("\nh =", h)
print("\ndelay_bar =", delay_bar)
print("\ngain_tap =", gain_tap)
print("\nQ =", Q)
print("\nM =", M)

fig, axes = plt.subplots(1, 2, figsize=(10,6))

axes[0].stem(h)
axes[0].set_xlabel("Índice da amostra n")
axes[0].set_ylabel("Amplitude")
axes[0].set_title("Resposta impulsiva discreta do canal")
axes[0].grid(True)

axes[1].stem(delay_bar, gain_tap)
axes[1].set_xlabel("Atraso (s)")
axes[1].set_ylabel("Amplitude (ganho)")
axes[1].set_title("Resposta do canal multipercurso")
axes[1].grid(True)

plt.tight_layout()
plt.show()

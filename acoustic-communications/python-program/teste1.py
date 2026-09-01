# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:15:08 2026

@author: Ricardo Alexandre
"""

import acousticChannel as ac

setup = {
    "Ts": 0.001,
    "paths": 6,
    "delayspread": 0.01}

gain = {
    "attenuation": 10}

delay = {
    "mean": 0.003}

doppler = {
    "type": "uniform",
    "velocity": 15}

h, delay_bar, gain_tap, Q, M = ac.acoustic_channel(setup, gain, delay, doppler)

print("h =", h)
print("delay_bar =", delay_bar)
print("gain_tap =", gain_tap)
print("Q =", Q)
print("M =", M)

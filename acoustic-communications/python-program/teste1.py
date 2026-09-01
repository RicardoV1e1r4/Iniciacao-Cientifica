# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:15:08 2026

@author: Ricardo Alexandre
"""

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

# -*- coding: utf-8 -*-
"""
Created on Sat Oct  10 17:52:29 2026

@author: Ricardo Alexandre
"""

import numpy as np
from scipy.stats import expon, rayleigh
from fractions import Fraction


def generate_ar1_gain_trajectory(num_samples, dt, f_d, sigma_gain=1.0):
    """Gera uma trajetória estocástica de ganho com autocorrelação temporal
    AR-1 baseada no espalhamento Doppler f_d.
    """
    rho = np.exp(-2 * np.pi * f_d * dt)
    w = (np.random.randn(num_samples) + 1j * np.random.randn(num_samples)) / np.sqrt(2)
    gain_complex = np.zeros(num_samples, dtype=complex)

    gain_complex[0] = w[0] * sigma_gain
    for k in range(1, num_samples):
        gain_complex[k] = rho * gain_complex[k - 1] + np.sqrt(1 - rho**2) * w[k] * sigma_gain

    return np.abs(gain_complex)


def method_b_taylor_coefficients(gain_minus, gain_zero, gain_plus, dt):
    """Calcula os coeficientes do polinômio de Taylor via diferenças
    finitas centradas em t = 0.
    """
    A_0 = gain_zero
    A_1 = (gain_plus - gain_minus) / (2 * dt)
    A_2 = (gain_plus - 2 * gain_zero + gain_minus) / (dt**2)

    return A_0, A_1, A_2


def acoustic_channel(setup, gain, delay, doppler):
    """
    Generate the impulse response of an underwater acoustic (UWA) channel.
    """

    # ==============================================================
    # 1. GENERATING THE DISTRIBUTION OF DELAYS
    # ==============================================================
    distribution_delay = expon(scale=delay["mean"])
    Dtau = distribution_delay.rvs(size=setup["paths"])
    Dtau_index = np.ceil(Dtau / setup["Ts"]).astype(int)
    delay_index = np.cumsum(Dtau_index)
    delay_values = (delay_index - 1) * setup["Ts"]

    # ==============================================================
    # 2. GENERATING THE GAINS (STOCHASTIC OR DETERMINISTIC/TAYLOR)
    # ==============================================================
    gain_mode = gain.get("mode", "stochastic")

    if gain_mode == "deterministic":
        # Parâmetros para geração determinística via Taylor
        f_d = doppler.get("f_d", 2.0)
        dt = gain.get("dt", 50 * setup["Ts"])
        T_bl = gain.get("T_bl", 0.025)

        # 1. Amostragem em t = -dt, 0, +dt para cada percurso
        gain_tap = np.zeros(setup["paths"])
        taylor_coeffs = []

        for p in range(setup["paths"]):
            dt_samples = generate_ar1_gain_trajectory(3, dt, f_d, sigma_gain=0.8)
            g_minus, g_zero, g_plus = dt_samples[0], dt_samples[1], dt_samples[2]

            A0, A1, A2 = method_b_taylor_coefficients(g_minus, g_zero, g_plus, dt)
            taylor_coeffs.append((A0, A1, A2))

            # Avalia o ganho médio ou instantâneo no bloco OFDM para o tap do canal
            t_vec = np.linspace(0, T_bl, max(1, int(T_bl / setup["Ts"])))
            g_taylor_vec = A0 + A1 * t_vec + 0.5 * A2 * (t_vec**2)
            gain_tap[p] = np.mean(g_taylor_vec)

    else:
        # Modo tradicional estocástico (Rayleigh)
        alpha = (
            np.log(10 ** (gain["attenuation"] / 10))
            / setup["delayspread"]
        )
        gain_variance = np.exp(-alpha * delay_values)
        rayleigh_scale = np.sqrt(gain_variance * 2 / (4 - np.pi))
        gain_tap = rayleigh.rvs(scale=rayleigh_scale)

    # Normalização dos ganhos
    gain_tap = gain_tap / np.linalg.norm(gain_tap)

    # ==============================================================
    # 3. DOPPLER EFFECT
    # ==============================================================
    c = 1500.0
    doppler_type = doppler["type"]

    if doppler_type == "none":
        a_max = 0.0
        Q, M = 1, 1
    elif doppler_type == "uniform":
        a_max = doppler["velocity"] / c
        fraction = Fraction(1 + a_max).limit_denominator()
        Q, M = fraction.numerator, fraction.denominator
    elif doppler_type == "non-uniform":
        a = np.asarray(doppler["velocity"]) / c
        a_max = np.max(a)
        fraction = Fraction(1 + a_max).limit_denominator()
        Q, M = fraction.numerator, fraction.denominator
    else:
        raise ValueError("Invalid Doppler type. Use 'none', 'uniform' or 'non-uniform'.")

    # ==============================================================
    # 4. CORRECTING THE DELAYS DUE TO DOPPLER
    # ==============================================================
    delay_index_bar = np.ceil(delay_index / (1 + a_max)).astype(int)
    delay_bar = (delay_index_bar - 1) * setup["Ts"]

    # ==============================================================
    # 5. GENERATING THE CHANNEL IMPULSE RESPONSE
    # ==============================================================
    h = np.zeros(np.max(delay_index_bar))
    h[delay_index_bar - 1] = gain_tap

    return h, delay_bar, gain_tap, Q, M
import numpy as np
import matplotlib.pyplot as plt


def ThorpsFormula(f):
	"""Calcula a atenuação acústica em dB/km para a frequência f em kHz."""
	f_squared = f**2
	return (
		0.11 * f_squared / (1 + f_squared)
		+ 44 * f_squared / (4100 + f_squared)
		+ 2.75e-4 * f_squared
		+ 0.003
	)


def conventionalPropagationLoss(f, d):
	"""Calcula a perda de propagação convencional em dB para a frequência f em kHz e distância d em km."""
	alpha = ThorpsFormula(f)
	return 20 * np.log10(d) + (alpha * d / 1000)


def transmissionLoss(f, d, L0 = 1, Ks = 1):
	"""Calcula a perda de transmissão em dB para a frequência f em kHz e distância d em km."""
	alpha = ThorpsFormula(f)
	TL = L0 * (d ** Ks) * (10**((alpha / 10000) * d))
	return TL


def distanceVariation(x_start, x_end, delta_x, y_start, y_end, delta_y):
	"""Gera um array de distâncias variando de d_start a d_end com incremento delta_d."""
	d_start = np.sqrt((x_start - x_end) ** 2 + (y_start - y_end) ** 2)
	d_end = np.sqrt((x_end + delta_x - x_start) ** 2 + (y_end + delta_y - y_start) ** 2)
	distance_variation = d_end - d_start
	return d_start, d_end, distance_variation

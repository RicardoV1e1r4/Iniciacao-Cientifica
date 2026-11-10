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
	return 20 * np.log10(d * 1000) + (alpha * d)


def transmissionLoss(f, d, L0 = 1, Ks = 1):
	"""Calcula a perda de transmissão em dB para a frequência f em kHz e distância d em km."""
	alpha = ThorpsFormula(f)
	TL = L0 * ((d * 1000) ** Ks) * (10**((alpha / 10) * d))
	return TL



def taylorAmplitudeVariation(A0, g_minus, g_plus, time_grid, order=1, g_zero=None):
	"""Calcula a variação temporal da amplitude A_p(t) de UM caminho por série de Taylor.

	Baseado na eq. (1.19) de Zhou & Wang (2014), seção 1.3.4. O tempo é centrado
	no meio do bloco (t em [-Tbl/2, +Tbl/2]) e a curva passa exatamente por
	A0*g_minus no início, A0*g_zero no centro e A0*g_plus no fim do bloco.

	Parâmetros:
		A0        : ganho de referência do percurso (escala linear, ex.: da grande escala).
		g_minus   : fator de ganho no início do bloco.
		g_plus    : fator de ganho no fim do bloco.
		time_grid : vetor de tempo dentro do bloco (ex.: np.linspace(0, Tbl, N)).
		order     : ordem da expansão (1 = linear, 2 = parabólica).
		g_zero    : fator de ganho no centro do bloco. Se None, usa a média de
		            g_minus e g_plus (na ordem 2 isso anula a curvatura).

	Retorno:
		A_t : vetor com a amplitude variante no tempo A_p(t).
	"""
	time_grid = np.asarray(time_grid, dtype=float)
	dt = time_grid[-1] - time_grid[0]
	if dt == 0:
		return np.full_like(time_grid, A0 * (g_zero if g_zero is not None else g_minus))

	# Tempo centrado no meio do bloco
	t_centered = time_grid - (time_grid[0] + dt / 2.0)

	g0 = g_zero if g_zero is not None else (g_minus + g_plus) / 2.0

	# Primeira derivada por diferença central (passo dt/2 para cada lado)
	A1 = (g_plus - g_minus) / dt

	if order == 1:
		A2 = 0.0
	elif order == 2:
		# Segunda derivada por diferença central com passo h = dt/2: (g+ - 2g0 + g-)/h^2
		A2 = (4.0 / dt**2) * (g_plus - 2.0 * g0 + g_minus)
	else:
		raise ValueError("A ordem deve ser 1 ou 2.")

	return A0 * (g0 + A1 * t_centered + 0.5 * A2 * t_centered**2)

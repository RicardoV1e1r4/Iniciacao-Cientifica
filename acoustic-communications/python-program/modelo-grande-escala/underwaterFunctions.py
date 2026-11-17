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


def conventionalPropagationLoss(f, d, Ks=2):
	"""Calcula a perda de propagação convencional em dB para a frequência f em kHz e distância d em km.

	Aceita escalares ou arrays (inclusive malhas de np.meshgrid). O cálculo é feito
	diretamente em dB, sem passar pela escala linear, portanto não há overflow.
	"""
	f = np.asarray(f, dtype=float)
	d = np.asarray(d, dtype=float)
	alpha = ThorpsFormula(f)
	return Ks * 10 * np.log10(d * 1000) + alpha * d


def transmissionLoss(f, d, L0=1, Ks=2, dB=True):
	"""Calcula a perda de transmissão para a frequência f em kHz e distância d em km.

	Modelo: TL = L0 * (d*1000)**Ks * 10**(alpha*d/10)

	O cálculo é feito em dB, evitando o overflow do valor linear em frequências e
	distâncias altas (ex.: 100 kHz a 100 km -> ~10^350):
		TL_dB = 10*log10(L0) + Ks*10*log10(d*1000) + alpha*d

	Parâmetros:
		dB : se True (padrão), retorna a perda em dB. Se False, retorna a perda em
		     escala linear (valores que excedem o limite do float viram inf).
	"""
	f = np.asarray(f, dtype=float)
	d = np.asarray(d, dtype=float)
	alpha = ThorpsFormula(f)
	TL_dB = 10 * np.log10(L0) + Ks * 10 * np.log10(d * 1000) + alpha * d
	if dB:
		return TL_dB
	with np.errstate(over='ignore'):
		return 10 ** (TL_dB / 10)



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


def chebyshevPoints(M, Tbl, t_start=0.0, lobatto=False):
	"""Retorna M instantes de amostragem de Chebyshev dentro do bloco [t_start, t_start + Tbl].

	Os pontos ficam mais concentrados perto das bordas do bloco, o que evita as
	oscilações do fenômeno de Runge em polinômios de ordem alta.

	Parâmetros:
		M       : número de pontos (para uma série de ordem N, use M = N + 1).
		Tbl     : duração do bloco em s.
		t_start : instante inicial do bloco em s.
		lobatto : se False, usa os pontos de Chebyshev clássicos
		          x_k = cos((2k - 1)*pi / (2M)), k = 1..M, que não incluem as bordas.
		          Se True, usa os pontos de Chebyshev-Lobatto x_k = cos(k*pi / (M - 1)),
		          k = 0..M-1, que incluem o início e o fim do bloco.

	Retorno:
		t_k : vetor com os M instantes, em ordem crescente.
	"""
	if M < 1:
		raise ValueError("M deve ser pelo menos 1.")
	if lobatto:
		if M < 2:
			raise ValueError("Os pontos de Chebyshev-Lobatto exigem M >= 2.")
		x_k = np.cos(np.arange(M) * np.pi / (M - 1))
	else:
		k = np.arange(1, M + 1)
		x_k = np.cos((2 * k - 1) * np.pi / (2 * M))
	# Converte de [-1, 1] para o bloco [t_start, t_start + Tbl]
	return np.sort(t_start + (x_k + 1) * Tbl / 2)


def taylorAmplitude(A0, gains, t_samples, time_grid, order=None, return_derivatives=False):
	"""Calcula a amplitude A_p(t) de UM caminho por série de Taylor de ordem qualquer.

	Generaliza a eq. (1.19) de Zhou & Wang (2014), seção 1.3.4:
		A_p(t) = sum_{n=0}^{N} A_p^(n) t^n / n!
	com t medido a partir do centro do bloco. As derivadas A_p^(n) são obtidas
	ajustando um polinômio de grau N aos ganhos medidos nos instantes t_samples.
	Com N + 1 ganhos, o polinômio passa exatamente por todos eles. Com mais
	ganhos, o ajuste é feito por mínimos quadrados.

	Internamente o tempo é normalizado para u em [-1, 1], o que evita a perda de
	precisão de potências altas como t^8.

	Parâmetros:
		A0                 : amplitude de referência do percurso (escala linear).
		gains              : ganhos relativos g_k = A(t_k) / A0 nos instantes t_samples.
		t_samples          : instantes, em s, em que os ganhos foram medidos.
		time_grid          : vetor de tempo do bloco onde a amplitude é calculada.
		                     O bloco vai de time_grid[0] a time_grid[-1].
		order              : ordem N da série. Se None, usa N = len(gains) - 1.
		return_derivatives : se True, também retorna as derivadas A_p^(n) no centro do bloco.

	Retorno:
		A_t         : vetor com a amplitude variante no tempo A_p(t).
		derivatives : (só se return_derivatives=True) vetor com A_p^(0), A_p^(1), ..., A_p^(N),
		              em unidades de amplitude por s^n.
	"""
	gains = np.asarray(gains, dtype=float)
	t_samples = np.asarray(t_samples, dtype=float)
	time_grid = np.asarray(time_grid, dtype=float)

	if gains.shape != t_samples.shape or gains.ndim != 1:
		raise ValueError("gains e t_samples devem ser vetores com o mesmo tamanho.")
	if order is None:
		order = len(gains) - 1
	if order < 0:
		raise ValueError("A ordem deve ser 0 ou maior.")
	if len(gains) < order + 1:
		raise ValueError(f"A ordem {order} exige pelo menos {order + 1} ganhos; foram dados {len(gains)}.")

	# Centro e meia-duração do bloco
	t_center = (time_grid[0] + time_grid[-1]) / 2.0
	half = (time_grid[-1] - time_grid[0]) / 2.0
	if half == 0:
		half = 1.0

	# Tempo normalizado: u = -1 no início do bloco e u = +1 no fim
	u_samples = (t_samples - t_center) / half
	u_grid = (time_grid - t_center) / half

	# Coeficientes c_n do polinômio sum c_n u^n (do grau 0 ao grau N)
	c = np.polynomial.polynomial.polyfit(u_samples, gains, order)
	A_t = A0 * np.polynomial.polynomial.polyval(u_grid, c)

	if not return_derivatives:
		return A_t

	# Derivadas no centro do bloco: A^(n) = A0 * n! * c_n / half^n
	n = np.arange(order + 1)
	factorials = np.array([np.prod(np.arange(1, k + 1), dtype=float) for k in n])
	derivatives = A0 * factorials * c / half**n
	return A_t, derivatives

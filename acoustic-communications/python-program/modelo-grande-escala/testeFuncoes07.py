import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf

# TESTE DA FUNÇÃO taylorAmplitude() (série de Taylor de ordem qualquer)
# Cenário: o receptor parte de (x_R0, y_R0) e o transmissor de (x_T0, y_T0).
# Cada um anda em linha reta com velocidade constante (v_R, v_T) na direção
# dada pelo seu ângulo (theta_R, theta_T), medido a partir do eixo x.
# A amplitude "exata" do caminho direto é relativa ao início do bloco:
#   A(t)/A(0) = 1/sqrt( (1 + dd/d0)^Ks * 10^(alpha*dd/10) ),  dd = d(t) - d0
# A série de Taylor de ordem N recebe os ganhos em M = N + 1 instantes do bloco
# (pontos de Chebyshev ou igualmente espaçados) e deve aproximar essa curva.


# Parâmetros estaticos
#-------------------------------------------------
f1      = 10        # Frequência da portadora em kHz
Ks      = 1.5       # Fator de espalhamento (1 = cilíndrico, 2 = esférico, 1.5 = prático)
v_R     = 10         # Velocidade do receptor em m/s
v_T     = 0         # Velocidade do transmissor em m/s
theta_R = 0         # Direção do receptor em graus (0 = +x, 90 = +y, 180 = -x)
theta_T = 0         # Direção do transmissor em graus
Tbl     = 20        # Duração do bloco em s (exagerada para a curvatura ficar visível)
N       = 1001      # Número de amostras de tempo dentro do bloco
ordem   = 4         # Ordem da série de Taylor
M       = ordem + 1 # Número de ganhos amostrados (uma série de ordem N exige N + 1 pontos)
ordem_max = 10      # Maior ordem testada na comparação de erros
#------------------------------------------------

# Componentes da velocidade do receptor, em m/s
v_Rx = v_R * np.cos(np.deg2rad(theta_R))
v_Ry = v_R * np.sin(np.deg2rad(theta_R))

# Componentes da velocidade do transmissor, em m/s
v_Tx = v_T * np.cos(np.deg2rad(theta_T))
v_Ty = v_T * np.sin(np.deg2rad(theta_T))

# Posição inicial do receptor, em km (as seguintes saem de v_R, theta_R e tempo)
x_R0 = 0
y_R0 = 0

# Posição inicial do transmissor, em km (as seguintes saem de v_T, theta_T e tempo)
x_T0 = 2
y_T0 = 0.4

# Parâmetros variáveis
#------------------------------------------------
tempo = np.linspace(0, Tbl, N)   # Vetor de tempo dentro do bloco, de 0 a Tbl
#------------------------------------------------


# AMPLITUDE EXATA (GRANDE ESCALA)
#------------------------------------------------
def distancia(t):
    """Distância transmissor-receptor em km no instante t (em s)."""
    x_R1 = x_R0 + v_Rx * t / 1000
    y_R1 = y_R0 + v_Ry * t / 1000
    x_T1 = x_T0 + v_Tx * t / 1000
    y_T1 = y_T0 + v_Ty * t / 1000
    return np.hypot(x_T1 - x_R1, y_T1 - y_R1)


def amplitude(t):
    """Amplitude do caminho direto, relativa ao início do bloco."""
    dd = distancia(t) - distancia(0)  # Variação da distância em relação ao início do bloco
    return 1 / np.sqrt(((1 + (dd / distancia(0))) ** Ks) * (10 ** (uwf.ThorpsFormula(f1) / 10 * dd)))


A_exata = amplitude(tempo)
print(f"Distância: início = {distancia(0):.3f} km | fim = {distancia(Tbl):.3f} km")


# SÉRIE DE TAYLOR DE ORDEM N
#------------------------------------------------
A0 = amplitude(Tbl / 2)                     # Amplitude de referência: centro do bloco

def aproximar(t_k, n):
    """Aproxima A(t) por uma série de ordem n usando os ganhos nos instantes t_k."""
    gains = amplitude(t_k) / A0             # Ganhos relativos g_k = A(t_k) / A0
    return uwf.taylorAmplitude(A0, gains, t_k, tempo, order=n, return_derivatives=True)

def erro_max(A_aprox):
    """Erro relativo máximo, em %, em relação à amplitude exata."""
    return 100 * np.max(np.abs(A_aprox - A_exata) / A_exata)

tk_cheb = uwf.chebyshevPoints(M, Tbl)       # Pontos de Chebyshev
tk_equi = np.linspace(0, Tbl, M)            # Pontos igualmente espaçados

A_cheb, derivadas = aproximar(tk_cheb, ordem)
A_equi, _ = aproximar(tk_equi, ordem)


# VERIFICAÇÕES NUMÉRICAS
#------------------------------------------------
# Com M = ordem + 1 pontos, o polinômio deve passar exatamente pelos ganhos usados
for t_k in (tk_cheb, tk_equi):
    A_nos = uwf.taylorAmplitude(A0, amplitude(t_k) / A0, t_k, t_k, order=ordem)
    assert np.allclose(A_nos, amplitude(t_k)), "A série não passa pelos pontos amostrados"

print(f"Pontos de Chebyshev (s): {np.round(tk_cheb, 2)}")
print(f"Derivadas no centro do bloco, A^(n) em amplitude/s^n:")
for n, An in enumerate(derivadas):
    print(f"    n = {n}: {An: .4e}")
print(f"Ordem {ordem}: erro máximo Chebyshev = {erro_max(A_cheb):.4f} % | "
      f"igualmente espaçados = {erro_max(A_equi):.4f} %")

# Erro máximo para várias ordens
ordens = np.arange(1, ordem_max + 1)
erros_cheb = [erro_max(aproximar(uwf.chebyshevPoints(n + 1, Tbl), n)[0]) for n in ordens]
erros_equi = [erro_max(aproximar(np.linspace(0, Tbl, n + 1), n)[0]) for n in ordens]


# GRÁFICOS
#------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]

ax1.plot(tempo, A_exata, 'k', linewidth=2, label='Exata')
ax1.plot(tempo, A_cheb, '--', label=f'Ordem {ordem}, Chebyshev')
ax1.plot(tempo, A_equi, '-.', label=f'Ordem {ordem}, igualmente espaçados')
ax1.plot(tk_cheb, amplitude(tk_cheb), 'o', label='Pontos de Chebyshev')
ax1.set_xlabel('Tempo (s)')
ax1.set_ylabel('A(t) / A(0)')
ax1.set_title('Amplitude de um caminho dentro do bloco')
ax1.legend()
ax1.grid(True, linestyle='--', linewidth=0.5)

ax2.plot(tempo, 100 * np.abs(A_cheb - A_exata) / A_exata, '--', label='Chebyshev')
ax2.plot(tempo, 100 * np.abs(A_equi - A_exata) / A_exata, '-.', label='Igualmente espaçados')
ax2.set_xlabel('Tempo (s)')
ax2.set_ylabel('Erro relativo (%)')
ax2.set_title(f'Erro da série de ordem {ordem}')
ax2.legend()
ax2.grid(True, linestyle='--', linewidth=0.5)

ax3.semilogy(ordens, erros_cheb, 'o-', label='Chebyshev')
ax3.semilogy(ordens, erros_equi, 's--', label='Igualmente espaçados')
ax3.set_xlabel('Ordem da série de Taylor')
ax3.set_ylabel('Erro relativo máximo (%)')
ax3.set_title('Erro em função da ordem')
ax3.legend()
ax3.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.show()

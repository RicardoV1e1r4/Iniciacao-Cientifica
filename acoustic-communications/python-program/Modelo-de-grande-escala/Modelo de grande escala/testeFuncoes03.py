import numpy as np
import matplotlib.pyplot as plt
import underwaterFunctions as uwf

# TESTE DA FUNÇÃO taylorAmplitudeVariation()
# Cenário: o receptor parte de (x_R0, y_R0) e o transmissor de (x_T0, y_T0).
# Cada um anda em linha reta com velocidade constante (v_R, v_T) na direção
# dada pelo seu ângulo (theta_R, theta_T), medido a partir do eixo x.
# A amplitude "exata" do caminho direto é A(t) = 1/sqrt(TL(d(t))), calculada
# com uwf.transmissionLoss(). A série de Taylor (ordens 1 e 2) só recebe os
# ganhos no início, no centro e no fim do bloco, e deve aproximar essa curva.


# Parâmetros estaticos
#-------------------------------------------------
f1 = 10         # Frequência da portadora em kHz
L0 = 1          # Perda de referência (escala linear)
Ks = 1.5        # Fator de espalhamento (1 = cilíndrico, 2 = esférico, 1.5 = prático)
v_R = 5           # Velocidade do receptor em m/s
v_T = 0           # Velocidade do transmissor em m/s
theta_R = 0       # Direção do receptor em graus (0 = +x, 90 = +y, 180 = -x)
theta_T = 0       # Direção do transmissor em graus
Tbl = 60        # Duração do bloco em s (exagerada para a curvatura ficar visível)
N = 1001        # Número de amostras de tempo dentro do bloco
#------------------------------------------------

# Componentes da velocidade do receptor, em m/s
v_Rx = v_R * np.cos(np.deg2rad(theta_R))
v_Ry = v_R * np.sin(np.deg2rad(theta_R))

# Componentes da velocidade do transmissor, em m/s
v_Tx = v_T * np.cos(np.deg2rad(theta_T))
v_Ty = v_T * np.sin(np.deg2rad(theta_T))

# Posição inicial do receptor, em km (a final sai de v_R, theta_R e Tbl)
x_R0 = 0
y_R0 = 0

# Posição inicial do transmissor, em km (a final sai de v_T, theta_T e Tbl)
x_T0 = 0.4
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
    """Amplitude do caminho direto: inverso da raiz da perda de transmissão."""
    dd = distancia(t) - distancia(0)  # Variação da distância em relação ao início do bloco
    # return 1 / np.sqrt(uwf.transmissionLoss(f1, distancia(t), L0, Ks))
    return 1 / np.sqrt(((1 + (dd / distancia(0))) ** Ks) * (10 ** (uwf.ThorpsFormula(f1) / 10 * dd)))


A_exata = amplitude(tempo)
print(f"Distância: início = {distancia(0):.3f} km | fim = {distancia(Tbl):.3f} km")


# PARÂMETROS DA SÉRIE DE TAYLOR
#------------------------------------------------
A0 = amplitude(Tbl / 2)                 # Amplitude de referência: centro do bloco
g_minus = amplitude(0) / A0             # Ganho relativo no início do bloco
g_zero = amplitude(Tbl / 2) / A0        # Ganho relativo no centro (vale 1)
g_plus = amplitude(Tbl) / A0            # Ganho relativo no fim do bloco

A_ordem1 = uwf.taylorAmplitudeVariation(A0, g_minus, g_plus, tempo, order=1)
A_ordem2 = uwf.taylorAmplitudeVariation(A0, g_minus, g_plus, tempo, order=2, g_zero=g_zero)


# VERIFICAÇÕES NUMÉRICAS
#------------------------------------------------
print(f"A0 = {A0:.4e} | g_minus = {g_minus:.4f} | g_zero = {g_zero:.4f} | g_plus = {g_plus:.4f}")

# A ordem 2 deve passar exatamente pelos três pontos usados no ajuste
for nome, k in [("início", 0), ("centro", N // 2), ("fim", -1)]:
    print(f"{nome:>6}: exata = {A_exata[k]:.4e} | ordem 2 = {A_ordem2[k]:.4e}")
    assert np.isclose(A_ordem2[k], A_exata[k]), f"Ordem 2 não passa pelo ponto do {nome}"

# A ordem 1 é uma reta: deve passar pelos pontos do início e do fim do bloco
assert np.isclose(A_ordem1[0], A_exata[0]), "Ordem 1 não passa pelo ponto do início"
assert np.isclose(A_ordem1[-1], A_exata[-1]), "Ordem 1 não passa pelo ponto do fim"

erro1 = 100 * np.abs(A_ordem1 - A_exata) / A_exata
erro2 = 100 * np.abs(A_ordem2 - A_exata) / A_exata
print(f"Erro relativo máximo: ordem 1 = {erro1.max():.3f} % | ordem 2 = {erro2.max():.3f} %")


# GRÁFICOS
#------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax1 = axes[0]
ax2 = axes[1]

ax1.plot(tempo, A_exata, 'k', linewidth=2, label='Exata (perda de transmissão)')
ax1.plot(tempo, A_ordem1, '--', label='Taylor ordem 1')
ax1.plot(tempo, A_ordem2, '-.', label='Taylor ordem 2')
ax1.plot([0, Tbl / 2, Tbl], A0 * np.array([g_minus, g_zero, g_plus]), 'ro', label='Pontos usados')
ax1.set_xlabel('Tempo (s)')
ax1.set_ylabel('Amplitude A_p(t)')
ax1.set_title('Amplitude de um caminho dentro do bloco')
ax1.legend()
ax1.grid(True, linestyle='--', linewidth=0.5)

ax2.plot(tempo, erro1, '--', label='Ordem 1')
ax2.plot(tempo, erro2, '-.', label='Ordem 2')
ax2.set_xlabel('Tempo (s)')
ax2.set_ylabel('Erro relativo (%)')
ax2.set_title('Erro da aproximação de Taylor')
ax2.legend()
ax2.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.show()

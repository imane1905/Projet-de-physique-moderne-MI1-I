from numpy import pi, exp, sqrt, real, imag, zeros, linspace, empty
import matplotlib.pyplot as plt

hbar = 1
m = 1

def GaussWP(k0, a_paquet, x, t):
    denominateur = m * a_paquet**2 + 2j * hbar * t
    val1 = sqrt(a_paquet) * (2*pi)**(-3/4) * sqrt(4 * pi * m / denominateur)
    exp1 = exp(1j * (k0 * x - hbar * k0**2 * t / (2*m)))
    exp2 = exp(-(x - hbar*k0*t/m)**2 / (a_paquet**2 + 2j*hbar*t/m))
    return val1*exp1*exp2

nx = 600
nt = 80000

x = linspace(-30, 30, nx)
t = linspace(0, 8, nt)

dx = x[1] - x[0]
dt = t[1] - t[0]

k0 = 5
a_paquet = 1
x0 = -20

seuil_probabilite = 1e-3   # seul changement par rapport à ton code d'origine

def simuler_et_mesurer(a_barriere, V0, X_cible, chercher_a_droite):
    V = zeros(nx)
    for i in range(nx):
        if 0 <= x[i] <= a_barriere:
            V[i] = V0

    psi = GaussWP(k0, a_paquet, x - x0, 0)

    tau = None
    for it in range(nt):
        somme_position = 0
        somme_densite = 0
        for i in range(nx):
            if (chercher_a_droite == False) or (x[i] > a_barriere):
                densite_i = abs(psi[i])**2
                somme_position = somme_position + x[i]*densite_i
                somme_densite = somme_densite + densite_i

        if somme_densite * dx > seuil_probabilite:
            position_moyenne = somme_position / somme_densite
            if position_moyenne >= X_cible:
                tau = t[it]
                break

        psi_nouveau = empty(nx, dtype=complex)
        psi_nouveau[0] = 0
        psi_nouveau[nx-1] = 0
        for i in range(1, nx - 1):
            d2psi = (psi[i+1] - 2*psi[i] + psi[i-1]) / dx**2
            psi_nouveau[i] = psi[i] + 1j * dt * ( (hbar/(2*m)) * d2psi - (V[i]/hbar) * psi[i] )
        psi = psi_nouveau

    return tau


tau0_num = simuler_et_mesurer(1, 0, 10, False)
print("tau0,num =", tau0_num)

liste_V0 = [6, 8, 10, 12, 13, 15, 18, 20]
liste_tau_t = []

for V0 in liste_V0:
    tau_t = simuler_et_mesurer(1, V0, 10, True)
    liste_tau_t.append(tau_t)
    print("V0 =", V0, " tau_t,num =", tau_t)

fig, ax = plt.subplots()
ax.plot(liste_V0, liste_tau_t, marker='o', label="tau_t,num (avec barriere)")
ax.axhline(y=tau0_num, color='gray', linestyle='--', label="tau0,num (cas libre)")
ax.set_xlabel("Hauteur de la barrière V0")
ax.set_ylabel("Temps de traversée")
ax.legend()
plt.show()
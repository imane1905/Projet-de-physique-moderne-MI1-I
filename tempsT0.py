import numpy as np
import matplotlib.pyplot as plt


hbar = 1
m = 1
# constantes 

Nx = 400
nt = 80000
# nombre de points en espace et nombre de pas de temps

x = np.linspace(-30, 30, Nx)
t = np.linspace(0, 8, nt)
# grille spatiale et grille temporelle

dx = x[1] - x[0]
dt = t[1] - t[0]
# pas d'espace et pas de temps

k0 = 5
x0 = -20
# vitesse initiale du paquet, position de départ

V0 = 0
a_barriere = 1
# ici V0=0, donc il n'y a en fait pas de vraie barrière, c'est le cas "libre"

V = np.zeros(Nx)
for i in range(Nx):
    if 0 <= x[i] <= a_barriere:
        V[i] = V0
# on construit V quand même comme dans les autres codes, mais comme V0=0
# ça ne change rien : c'est juste pour garder la même structure de code

# Initialisation
psi = np.exp(-((x - x0)**2) / 4) * np.exp(1j * k0 * x)
# paquet d'ondes gaussien initial, centré en x0, avec une vitesse k0

for it in range(nt - 1):
    for ix in range(1, Nx - 1):
        # dérivée seconde
        derivee_seconde = (psi[ix+1] - 2*psi[ix] + psi[ix-1]) / (dx**2)

        # schrodinger
        energie_cinetique = - (hbar**2 / (2 * m)) * derivee_seconde
        energie_potentielle = V[ix] * psi[ix]
        dpsi_dt = (energie_cinetique + energie_potentielle) / (1j * hbar)

        # Euler
        psi[ix] = psi[ix] + dt * dpsi_dt

    if it == 0 or it == nt//3 or it == 2*nt//3 or it == nt-2:
        plt.plot(x, np.abs(psi)**2, label=f"t={it*dt:.2f}")
        # on trace la densité de probabilité à quelques instants choisis
        # (début, un tiers, deux tiers, fin)

plt.title("Évolution du paquet d'ondes (Cas libre V0=0)")
plt.legend()
plt.show()
# titre, légende, affichage

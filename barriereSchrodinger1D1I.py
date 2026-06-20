import numpy as np
import matplotlib.pyplot as plt


hbar = 1
m = 1

Nx = 400
nt = 80000 

x = np.linspace(-30, 30, Nx)
t = np.linspace(0, 8, nt)
dx = x[1] - x[0]
dt = t[1] - t[0]

k0 = 5
x0 = -20

# Paramètres barrière 
V0 = 13.0
a_barriere = 1.0

V = np.zeros(Nx)
for i in range(Nx):
    if 0 <= x[i] <= a_barriere:
        V[i] = V0

# Initialisation (Paquet complexe)
psi = np.exp(-((x - x0)**2) / 4) * np.exp(1j * k0 * x)


fig, ax = plt.subplots()

for it in range(nt - 1):
    for ix in range(1, Nx - 1):
        # 1. Calcul de dérivée seconde
        derivee_seconde = (psi[ix+1] - 2*psi[ix] + psi[ix-1]) / (dx**2)
        
        # 2. Schrödinger
        energie_cinetique = - (hbar**2 / (2 * m)) * derivee_seconde
        energie_potentielle = V[ix] * psi[ix]
        dpsi_dt = (energie_cinetique + energie_potentielle) / (1j * hbar)
        
        # 3. Mise à jour Euler
        psi[ix] = psi[ix] + dt * dpsi_dt

    # Affichage aux instants clés
    if it == 0 or it == nt//4 or it == nt//2 or it == nt-2:
        ax.plot(x, np.abs(psi)**2, label=f"t={it*dt:.2f}")

# Dessin de la barrière
ax.axvspan(0, a_barriere, color='gray', alpha=0.3, label="Barrière V0=13")
ax.set_title(f"Effet Tunnel : a = {a_barriere}, V0 = {V0}")
ax.legend()
plt.show()
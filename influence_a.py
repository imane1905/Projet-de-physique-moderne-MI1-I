from numpy import pi, exp, sqrt, zeros, linspace, roll
import numpy as np
import matplotlib.pyplot as plt

hbar = 1
m = 1
k0 = 5
a_paquet = 1
x0 = -10
V0 = 13.0

nx = 400
nt = 80000
x = linspace(-30, 30, nx)
t = linspace(0, 8, nt)
dx = x[1] - x[0]
dt = t[1] - t[0]

seuil = 0.02  # 2% de la hauteur du paquet

def GaussWP(k0, a_paquet, x, t):
    denominateur = m * a_paquet**2 + 2j * hbar * t
    val1 = sqrt(a_paquet) * (2*pi)**(-3/4) * sqrt(4 * pi * m / denominateur)
    exp1 = exp(1j * (k0 * x - hbar * k0**2 * t / (2*m)))
    exp2 = exp(-(x - hbar*k0*t/m)**2 / (a_paquet**2 + 2j*hbar*t/m))
    return val1*exp1*exp2

def evoluer(psi, V):
    derivee_seconde = (roll(psi, -1) - 2*psi + roll(psi, 1)) / dx**2
    dpsi_dt = (-(hbar**2/(2*m))*derivee_seconde + V*psi) / (1j*hbar)
    psi_nouveau = psi + dt * dpsi_dt
    psi_nouveau[0] = 0
    psi_nouveau[-1] = 0
    return psi_nouveau

# Vitesse de groupe mesurée sur une simulation libre de référence
def mesurer_vg():
    psi = GaussWP(k0, a_paquet, x - x0, 0)
    V = zeros(nx)
    pos1 = pos2 = t1 = t2 = None
    for it in range(nt):
        if pos1 is None and t[it] >= 0.2:
            pos1, t1 = x[np.argmax(np.abs(psi)**2)], t[it]
        if pos2 is None and t[it] >= 1.0:
            pos2, t2 = x[np.argmax(np.abs(psi)**2)], t[it]
            break
        psi = evoluer(psi, V)
    return (pos2 - pos1) / (t2 - t1)

v_g = mesurer_vg()
t_entree = abs(x0) / v_g

def mesurer_tau0(a_barriere):
    psi = GaussWP(k0, a_paquet, x - x0, 0)
    V = zeros(nx)
    t_in = t_out = None
    for it in range(nt):
        position_pic = x[np.argmax(np.abs(psi)**2)]
        if t_in is None and position_pic >= 0:
            t_in = t[it]
        elif t_in is not None and t_out is None and position_pic >= a_barriere:
            t_out = t[it]
            break
        psi = evoluer(psi, V)
    return None if t_in is None or t_out is None else t_out - t_in

def mesurer_tau_t(a_barriere):
    V = zeros(nx)
    V[(x >= 0) & (x <= a_barriere)] = V0
    psi = GaussWP(k0, a_paquet, x - x0, 0)
    for it in range(nt):
        if t[it] >= t_entree:
            densite = np.abs(psi)**2
            masque = x > a_barriere
            if np.any(masque) and np.max(densite[masque]) > seuil * np.max(densite):
                return t[it] - t_entree
        psi = evoluer(psi, V)
    return None

liste_a = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
tau0_liste = [mesurer_tau0(a) for a in liste_a]
tau_t_liste = [mesurer_tau_t(a) for a in liste_a]

for a, t0, tt in zip(liste_a, tau0_liste, tau_t_liste):
    print(f"a = {a} -> tau0,num = {t0:.3f}, tau_t,num = {tt:.3f}")

fig, ax = plt.subplots()
ax.plot(liste_a, tau0_liste, 'go-', linewidth=2, label="τ0,num (libre)")
ax.plot(liste_a, tau_t_liste, 'bo-', linewidth=2, label="τt,num (barrière V0=13)")
ax.set_title("Influence de la largeur a sur τ0,num et τt,num")
ax.set_xlabel("Largeur de la barrière a")
ax.set_ylabel("Temps de traversée")
ax.legend()
plt.show()
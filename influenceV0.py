from numpy import pi, exp, sqrt, real, imag, zeros, linspace, empty
import matplotlib.pyplot as plt
# numpy pour les calculs (avec empty en plus cette fois), matplotlib pour le graphique

hbar = 1
m = 1
# constantes    

def GaussWP(k0, a_paquet, x, t):
    denominateur = m * a_paquet**2 + 2j * hbar * t
    val1 = sqrt(a_paquet) * (2*pi)**(-3/4) * sqrt(4 * pi * m / denominateur)
    exp1 = exp(1j * (k0 * x - hbar * k0**2 * t / (2*m)))
    exp2 = exp(-(x - hbar*k0*t/m)**2 / (a_paquet**2 + 2j*hbar*t/m))
    return val1*exp1*exp2
    # même fonction que dans le code du paquet d'ondes gaussien, pour initialiser Ψ

nx = 600
nt = 80000
x = linspace(-30, 30, nx)
t = linspace(0, 8, nt)
dx = x[1] - x[0]
dt = t[1] - t[0]
# grilles en espace et en temps comme avant

k0 = 5
a_paquet = 1
x0 = -20
# vitesse initiale, largeur, et position de départ du paquet

seuil_probabilite = 1e-3 
# seuil minimum de probabilité pour considérer la mesure comme fiable


def simuler_et_mesurer(a_barriere, V0, X_cible, chercher_a_droite):
    V = zeros(nx)
    for i in range(nx):
        if 0 <= x[i] <= a_barriere:
            V[i] = V0
    # barriere potentiel

    psi = GaussWP(k0, a_paquet, x - x0, 0)
    # on initialise le paquet d'ondes centré en x0

    tau = None
    # tau va stocker le temps qu'on cherche (temps de traversée), None si pas encore trouvé

    for it in range(nt):
        somme_position = 0
        somme_densite = 0
        for i in range(nx):
            if (chercher_a_droite == False) or (x[i] > a_barriere):
                densite_i = abs(psi[i])**2
                somme_position = somme_position + x[i]*densite_i
                somme_densite = somme_densite + densite_i
                # on calcule la position moyenne du paquet
                # si chercher_a_droite est True, on ne regarde que la partie qui a passé la barrière

        if somme_densite * dx > seuil_probabilite:
            position_moyenne = somme_position / somme_densite
            if position_moyenne >= X_cible:
                tau = t[it]
                break
                # si assez de probabilité ET paquet a atteint la position cible X_cible,
                # alors on note le temps actuel comme tau, et on arrête la boucle

        psi_nouveau = empty(nx, dtype=complex)
        psi_nouveau[0] = 0
        psi_nouveau[nx-1] = 0
        for i in range(1, nx - 1):
            d2psi = (psi[i+1] - 2*psi[i] + psi[i-1]) / dx**2
            psi_nouveau[i] = psi[i] + 1j * dt * ( (hbar/(2*m)) * d2psi - (V[i]/hbar) * psi[i] )
        psi = psi_nouveau
        # sinon, on continue à faire évoluer Ψ dans le temps avec Schrödinger,
            

    return tau
    # la fonction renvoie le temps de traversée mesuré (ou None si jamais atteint)

tau0_num = simuler_et_mesurer(1, 0, 10, False)
print("tau0,num =", tau0_num)
# sans barrière (V0=0), pour avoir un temps de traversée "libre"

liste_V0 = [6, 8, 10, 12, 13, 15, 18, 20]
liste_tau_t = []

for V0 in liste_V0:
    tau_t = simuler_et_mesurer(1, V0, 10, True)
    liste_tau_t.append(tau_t)
    print("V0 =", V0, " tau_t,num =", tau_t)
    # on teste plusieurs hauteurs de barrière, et pour chacune on mesure le temps que met le paquet pour atteindre X_cible
    

fig, ax = plt.subplots()
ax.plot(liste_V0, liste_tau_t, marker='o', label="tau_t,num (avec barriere)")
ax.axhline(y=tau0_num, color='gray', linestyle='--', label="tau0,num (cas libre)")
ax.set_xlabel("Hauteur de la barrière V0")
ax.set_ylabel("Temps de traversée")
ax.legend()
plt.show()
# on trace le temps de traversée en fonction de V0, avec une ligne de référence

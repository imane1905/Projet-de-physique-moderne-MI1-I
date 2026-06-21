import numpy as np
import matplotlib.pyplot as plt


hbar = 1
m = 1
# constantes 

V0 = 13.0
a = 1.0

E = np.linspace(0.01, 25, 1000)
# on fait varier l'énergie E entre 0.01 et 25, sur 1000 points


T = np.zeros_like(E)
# on prépare un tableau vide pour stocker le coefficient de transmission, même taille que E

for i in range(len(E)):
    if E[i] < V0:
        kappa = np.sqrt(2*m*(V0 - E[i])) / hbar
        T[i] = 1 / (1 + (V0**2 * np.sinh(kappa*a)**2) / (4*E[i]*(V0 - E[i])))
# cas où l'énergie est inférieure à la barrière régime tunnel (formule analytique)
    else:
        k2 = np.sqrt(2*m*(E[i] - V0)) / hbar
        T[i] = 1 / (1 + (V0**2 * np.sin(k2*a)**2) / (4*E[i]*(E[i] - V0)))
# cas où l'énergie dépasse la barrière la particule passe 

fig, ax = plt.subplots()
ax.plot(E, T, linewidth=2)
# on trace T en fonction de E
ax.axvline(x=V0, color='gray', linestyle='--', label=f"V0 = {V0}")
ax.scatter([12.5], [0.097], color='red', zorder=5, label="Notre paquet (E=12.5)")
# on place un point rouge pour montrer où se situe notre paquet d'ondes du code précédent
ax.set_xlabel("Énergie E")
ax.set_ylabel("Coefficient de transmission T")
ax.set_title("Coefficient de transmission en fonction de l'énergie")
ax.legend()
plt.show()
# noms des axes, titre, légende, affichage

import numpy as np
import matplotlib.pyplot as plt

hbar = 1
m = 1
V0 = 13.0
a = 1.0

E = np.linspace(0.01, 25, 1000)

T = np.zeros_like(E)

for i in range(len(E)):
    if E[i] < V0:
        kappa = np.sqrt(2*m*(V0 - E[i])) / hbar
        T[i] = 1 / (1 + (V0**2 * np.sinh(kappa*a)**2) / (4*E[i]*(V0 - E[i])))
    else:
        k2 = np.sqrt(2*m*(E[i] - V0)) / hbar
        T[i] = 1 / (1 + (V0**2 * np.sin(k2*a)**2) / (4*E[i]*(E[i] - V0)))

fig, ax = plt.subplots()
ax.plot(E, T, linewidth=2)
ax.axvline(x=V0, color='gray', linestyle='--', label=f"V0 = {V0}")
ax.scatter([12.5], [0.097], color='red', zorder=5, label="Notre paquet (E=12.5)")
ax.set_xlabel("Énergie E")
ax.set_ylabel("Coefficient de transmission T")
ax.set_title("Coefficient de transmission en fonction de l'énergie")
ax.legend()
plt.show()
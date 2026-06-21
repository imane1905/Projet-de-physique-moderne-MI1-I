from numpy import pi, exp, sqrt, real, imag, zeros, linspace
import matplotlib.pyplot as plt

hbar=1
m=1
# constantes

def GaussWP(k0, a, x, t):
    denominateur = m * a**2 + 2j * hbar * t
    val1 = sqrt(a) * (2*pi)**(-3/4) * sqrt(4 * pi * m / denominateur)
    # facteur de normalisation du paquet d'ondes
    exp1 = exp(1j * (k0 * x - hbar * k0**2 * t / (2*m)))
    # partie oscillante
    exp2 = exp(-(x - hbar*k0*t/m)**2 / (a**2 + 2j*hbar*t/m))
    # enveloppe gaussienne qui se déplace et s'étale avec le temps
    return val1*exp1*exp2
    # formule finale = multiplication des trois termes

k0 = 1
a = 1
t = 2

x = linspace(-10, 10, 300)
psi = GaussWP(k0, a, x, t)
fig, ax = plt.subplots()
ax.plot(x, real(psi), label="Re[Ψ(x,t)]")
ax.plot(x, imag(psi), label="Im[Ψ(x,t)]")
# on trace partie réelle et imaginaire
ax.set_title(f"Paquet d'ondes gaussien à t={t}")
ax.set_xlabel("x")
ax.set_ylabel("Ψ(x,t)")
ax.legend()
plt.show()





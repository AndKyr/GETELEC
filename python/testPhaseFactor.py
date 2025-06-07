import numpy as np

Ntest = 8
Nv = 16
v = np.linspace(0.,1.999 * np.pi, Nv)

u = np.random.rand(Ntest)
uPrime = np.random.rand(Ntest)
vPrime = np.exp(-2*u)

c1 = np.random.rand(Ntest) + 1j * np.random.rand(Ntest)
c2 = np.random.rand(Ntest) + 1j * np.random.rand(Ntest)

for i in range(Ntest):
    print("\n test ", i)
    for j in range(Nv):
        psi = np.exp(u[i]) * (np.cos(v[j]) + 1j * np.sin(v[j]))
        psiPrime = psi * (uPrime[i] + 1j * vPrime[i])

        PhiMatrix = np.array([[np.real(psi), np.imag(psi)],[np.real(psiPrime), np.imag(psiPrime)]])
        c = np.array([c1[i], c2[i]])
        result = PhiMatrix @ c
        print(np.angle(result[0]), np.angle(result[1]))
        



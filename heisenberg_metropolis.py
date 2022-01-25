"""
    3D heisenberg model simultion using Metropolis sampling
    João Inácio, jan 2022
"""

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

# System parameters
L = 4
N_atm = L**3
J = 1
nn = 6
nn_table = np.loadtxt("NN_tables/NN_3D_L" + str(L) + ".txt")

# Metropolis parameters
temperatures = np.linspace(0.01, 3, 10)
steps = int(5e4)
measure = int(4e4)
n_measure = steps - measure
print(f"temperatures = {temperatures}")

M = np.zeros(len(temperatures))
E = np.zeros(len(temperatures))

# Initialization of configuration and RNG
rng = npr.default_rng(npr.MT19937())
rng1 = npr.default_rng(npr.MT19937(seed=1))
rng2 = npr.default_rng(npr.MT19937(seed=2))

spin_lattice = dict()
for i in range(N_atm):
    x1 = rng1.uniform(-1, 1)
    x2 = rng2.uniform(-1, 1)
    while x1**2 + x2**2 >= 1:
        x1 = rng1.uniform(-1, 1)
        x2 = rng2.uniform(-1, 1)
    
    spin_lattice[i] = (2 * x1 * np.sqrt(1 - x1**2 - x2**2), 
                       2 * x2 * np.sqrt(1 - x1**2 - x2**2), 
                       1 - 2 * (x1**2 + x2**2))

# Metropolis Sampling
for t_idx, T in enumerate(temperatures):
    M_tmp = np.zeros(n_measure)
    E_tmp = np.zeros(n_measure)
    measure_idx = 0
    
    for t in range(steps):
        x1 = rng1.uniform(-1, 1)
        x2 = rng2.uniform(-1, 1)
        while x1**2 + x2**2 >= 1:
            x1 = rng1.uniform(-1, 1)
            x2 = rng2.uniform(-1, 1)
            
        new_spin = (2 * x1 * np.sqrt(1 - x1**2 - x2**2), 
                    2 * x2 * np.sqrt(1 - x1**2 - x2**2), 
                    1 - 2 * (x1**2 + x2**2))
        
        site = rng.integers(0, N_atm, endpoint=False)
        
        sum_NNx, sum_NNy, sum_NNz = 0, 0, 0
        for a in range(nn):
            sum_NNx += spin_lattice[nn_table[site, a]][0]
            sum_NNx += spin_lattice[nn_table[site, a]][1]
            sum_NNx += spin_lattice[nn_table[site, a]][2]
        
        delta_E = - J * ((new_spin[0] - spin_lattice[site][0]) * sum_NNx + 
                        (new_spin[1] - spin_lattice[site][1]) * sum_NNy + 
                        (new_spin[2] - spin_lattice[site][2]) * sum_NNz)
        
        if delta_E <= 0 or rng.random() <= np.exp(- delta_E / T):
            spin_lattice[site] = new_spin
        
        if t > measure:
            # Magnetization
            Mx, My, Mz = 0, 0, 0
            for i in range(N_atm):
                Mx += spin_lattice[i][0]
                My += spin_lattice[i][1]
                Mz += spin_lattice[i][2]
            M_tmp[measure_idx] = np.sqrt(Mx**2 + My**2 + Mz**2) / N_atm
            
            # Energy
            for i in range(N_atm):
                for a in range(nn):
                    E_tmp[measure_idx] += np.dot(spin_lattice[site], spin_lattice[nn_table[site, a]])
            E_tmp[measure_idx] = - J * E_tmp[measure_idx] / (2 * N_atm)
            
            measure_idx += 1

    M[t_idx] = np.mean(M_tmp)
    E[t_idx] = np.mean(E_tmp)
    
    print(f"T = {t_idx + 1}/{len(temperatures)}", end="\r")

print(f"M = {M}")
print(f"E = {E}")

# plt.figure(1)
# plt.plot(temperatures, M, 'o-')
# plt.xlabel(r"$T$")
# plt.ylabel(r"$M$")
# plt.title(r"$\langle M(T) \rangle$")

# plt.figure(2)
# plt.plot(temperatures, E, 'o-')
# plt.xlabel(r"$T$")
# plt.ylabel(r"$E$")
# plt.title(r"$\langle E(T) \rangle$")

# plt.show()



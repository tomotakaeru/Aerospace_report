
#%%
import numpy as np
import pprint

beta = np.arange(0, 180, 0.1) / 180 * np.pi
Vin = np.arange(0, 30, 0.1)

mu_x = 3.986e+15  #[m3/s2]
Vx = 35.60  #[km/s]
Vtx = 38.61  #[km/s]
rx = 5000  #[km]

e_i = 0.5 / 180 * np.pi
e_dV = 21.95
solve = []

for b in beta:
    for v in Vin:
        Vrev = np.sqrt(Vx ** 2 + v ** 2 - 2 * Vx * v * np.cos(b))
        a = 2 * np.arcsin(mu_x / (mu_x + (Vrev * 1000)** 2 * (rx * 1000)))
        a_dash = a - np.arcsin(v / Vrev * np.sin(b))
        i = np.arctan2(Vrev * np.sin(a_dash), Vx - Vrev * np.cos(a_dash))
        if abs(i - np.pi / 2) < e_i and np.sqrt(Vtx ** 2 + v ** 2 - 2 * Vtx * v * np.cos(b)) < e_dV:
            # e = abs(i - np.pi / 2)
            dV = np.sqrt(Vtx ** 2 + v ** 2 - 2 * Vtx * v * np.cos(b))
            solve.append(["{:>6.2f}".format(b / np.pi * 180), "{:>6.2f}".format(v), "{:>6.2f}".format(dV), "{:>6.2f}".format(i / np.pi * 180), "{:>6.2f}".format(a / np.pi * 180)])

print("[ beta, Vin, dV, i, alpha ]")
pprint.pprint(solve, width=70, compact=True)

#%%

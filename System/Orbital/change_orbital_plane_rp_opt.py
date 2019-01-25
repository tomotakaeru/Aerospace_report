
#%%
import numpy as np
import pprint

beta = np.arange(0, 180, 0.1) / 180 * np.pi
Vin = np.arange(0, 40, 0.1)

mu_x = 3.986e+15  #[m3/s2]
mu_s = 1.328e+20  #[m3/s2]
Vx = 35.60  #[km/s]
Vtx = 38.61  #[km/s]
rx = 5000  #[km]
Rx = 0.7 * 1.496e+11  #[m]

e_rp = 700000  #[km]
e_dV = 35  #[km/s]
solve = []

for b in beta:
    for v in Vin:
        Vrev = np.sqrt(Vx ** 2 + v ** 2 - 2 * Vx * v * np.cos(b))
        a = 2 * np.arcsin(mu_x / (mu_x + (Vrev * 1000)** 2 * (rx * 1000)))
        # a_dash = -np.pi + a + np.arcsin(v / Vrev * np.sin(b))
        a_dash = np.pi + a - np.arcsin(v / Vrev * np.sin(b))
        # i = np.arctan2(Vrev * np.sin(a_dash), Vx + Vrev * np.cos(a_dash))
        Vout = np.sqrt((Vrev * np.sin(a_dash)) ** 2 + (Vx + Vrev * np.cos(a_dash)) ** 2)
        rp = (-Rx - 1 / ((Vout * 1000)** 2 / (2 * mu_s) - 1 / Rx)) / 1000

        if abs(rp) < e_rp:
            dV = np.sqrt(Vtx ** 2 + v ** 2 - 2 * Vtx * v * np.cos(b))
            if dV < e_dV:
                # e_rp = rp
                e_dV = dV
                solve.append(["{:>6.2f}".format(b / np.pi * 180), "{:>6.2f}".format(dV), "{:>6.2f}".format(a / np.pi * 180), "{:>6.2f}".format(Vout), "{:>6.0f}".format(rp)])

print(" [   beta ,      dV ,   alpha ,    Vout ,  rp[km] ]")
pprint.pprint(solve, width=60, compact=True)
print(" [   beta ,      dV ,   alpha ,    Vout ,  rp[km] ]")

#%%

#3dos_eq_input_roll.pyの初期速度を200m/sにしてロール角も0で一定とした
#変数endにより，時間範囲を変えられるようにした

#指定されたロール角入力時の速度・位置・高度の推移

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m = 7500 #kg
S = 27.9 #m2
CLa = 4.30
CD0 = 0.0548
K = 3.02
rho = 0.736 #kg/m3
U = 200 #m/s
h = 5000 #m
g = 9.807 #m/s2

alpha = 0.0507 #rad
T = 20800 #N

end = 300 #何秒後まで見るか
n = 3001 #要素数
t_list =[i*end/(n-1) for i in range(0, n)]
# phi_list = [0 for i in t_list[:math.ceil((n-1)/80*20)+1]] + [np.pi/6*np.sin(np.pi/20*(i-20))
#  for i in t_list[math.ceil((n-1)/80*20)+1:math.ceil((n-1)/80*60)+1]] + [0 for i in t_list[math.ceil((n-1)/80*60)+1:]]
phi_list=[0 for i in range(0,n)]

U_list, gamma_list, psi_list = [U for i in range(n)], [0 for i in range(n)], [0 for i in range(n)]

def U_eq(i):
    return U_list[i] + ((- rho/2*pow(U_list[i],2)*S*(CD0+K*pow(alpha,2)) + T*np.cos(alpha))/m - g*np.sin(gamma_list[i])) * end/(n-1)
def gamma_eq(i):
    return gamma_list[i] + ((rho/2*pow(U_list[i],2)*S*CLa*alpha + T*np.sin(alpha))/m*np.cos(phi_list[i]) - g*np.cos(gamma_list[i]))/U_list[i] * end/(n-1)
def psi_eq(i):
    return psi_list[i] + ((rho/2*pow(U_list[i],2)*S*CLa*alpha + T*np.sin(alpha))/m*np.sin(phi_list[i]) / np.cos(gamma_list[i]))/U_list[i] * end/(n-1)

for i in range(n-1):
    U_list[i+1] = U_eq(i)
    gamma_list[i+1] = gamma_eq(i)
    psi_list[i+1] = psi_eq(i)


x_list, y_list, h_list = [0 for i in range(n)], [0 for i in range(n)], [h for i in range(n)]

def x_eq(i):
    return x_list[i] + U_list[i]*np.cos(gamma_list[i])*np.cos(psi_list[i]) * end/(n-1)
def y_eq(i):
    return y_list[i] + U_list[i]*np.cos(gamma_list[i])*np.sin(psi_list[i]) * end/(n-1)
def h_eq(i):
    return h_list[i] + U_list[i]*np.sin(gamma_list[i]) * end/(n-1)

for i in range(n-1):
    x_list[i+1] = x_eq(i)
    y_list[i+1] = y_eq(i)
    h_list[i+1] = h_eq(i)


"""
plt.xlabel("t [sec]")
plt.ylabel("U [m/s]")
plt.ylim(160,210)
plt.plot(t_list, U_list)
plt.show()
"""
plt.xlabel("t [sec]")
plt.ylabel("h [m]")
plt.ylim(4950,5650)
plt.plot(t_list, h_list)
plt.show()
"""
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.plot(x_list, y_list)
plt.show()
"""
"""
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
ax.scatter(x_list,y_list,h_list)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_zlabel("z [m]")
plt.show()
"""

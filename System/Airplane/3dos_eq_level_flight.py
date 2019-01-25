#以下の諸元の航空機が水平定常直線飛行を行うときの迎角と推力を求める

import numpy as np

m = 7500 #kg
S = 27.9 #m2
CLa = 4.30
CD0 = 0.0548
K = 3.02
rho = 0.736 #kg/m3
U = 180 #m/s
g = 9.807 #m/s2

def vertical_eq(alpha, T):
    return rho/2*pow(U,2)*S*CLa*alpha + T*np.sin(alpha)-m*g 
def horizontal_eq(alpha, T):
    return T*np.cos(alpha) - rho/2*pow(U,2)*S*(CD0+K*pow(alpha,2))

alpha_list = np.linspace(0, 0.2, 2001) #rad
T_list = np.linspace(10000, 80000, 701) #N

vertical_result=1.0
horizontal_result=1.0

for i in alpha_list:
    for j in T_list:
        vertical_result = vertical_eq(i, j)
        horizontal_result = horizontal_eq(i, j)
        if abs(vertical_result)<50 and abs(horizontal_result)<50:
            print("alpha:{}, T:{}".format(i,j))

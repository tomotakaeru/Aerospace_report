"""
ジュウコフスキー翼の作成
指定した迎角の翼の，内部の点と枠の点の座標を得る
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
from decimal import Decimal, ROUND_HALF_UP


# z平面の円から座標変換により翼型を生成
#初期条件
a = 0.6
gamma = 5.0
v0 = 1.0
rad = np.pi / 180.0
alpha = 5.0 * rad
beta = 20.0 * rad
c = 0.5

#計算範囲及び精度
n = 181
m = 181
dr = 0.015
dtheta = 2.0 * rad
x = np.zeros((n, n))
y = np.zeros((n, n))
X = np.zeros((n, n))
Y = np.zeros((n, n))
Xze = np.zeros((n, n))
Yze = np.zeros((n, n))
stream = np.zeros((n, n))
poten = np.zeros((n, n))
gamma = 4.0 * np.pi * v0 * a * np.sin(alpha + beta)  #Kutta条件
Zc = a * np.exp(1j * (np.pi - beta)) + c
cgamma = 1j * gamma / (2.0 * np.pi)

def potential(s, t):
    z = s + 1j * t
    f = v0 * (z + a ** 2 / z) + cgamma * np.log(z)
    return f
def circle_fc(s, t):
    z = s + 1j * t
    Z = z * np.exp(1j * alpha) + Zc
    return Z
def Zeta(s, t):
    z = s + 1j * t
    Z = z * np.exp(1j * alpha) + Zc
    Zeta = Z + c ** 2 / Z
    return Zeta
#z座標系からζ座標系の圧力分布に変換する関数
def pres(s, t):
    z = s + 1j * t
    Z = z * np.exp(1j * alpha) + Zc
    f = np.exp(-1j * alpha) * v0 * (1 - a ** 2 / z ** 2 + 1j * gamma / (v0 * 2 * np.pi * z)) / (1 - c ** 2 / Z ** 2)
    Cp = 1 - np.abs(f)** 2 / v0 ** 2
    return Cp

#変換
for i in range(n):
    for k in range(m):
        #z平面上で同心円状に計算点を配置
        theta0 = 0.0
        theta0 = -beta
        rrr = a + dr * i
        theta = theta0 + dtheta * k
        x[i, k] = rrr * np.cos(theta)
        y[i, k] = rrr * np.sin(theta)
        poten[i, k] = np.real(potential(x[i, k], y[i, k]))
        stream[i, k] = np.imag(potential(x[i, k], y[i, k]))
        Z = circle_fc(x[i, k], y[i, k])
        Xze[i, k] = np.real(Z + c ** 2.0 / Z)
        Yze[i, k] = np.imag(Z + c ** 2.0 / Z)
Wing = Zeta(x[0,:], y[0,:])  #翼表面部分にi=0が対応


# 迎角の指定
degree = 20
#迎角を0にする
Wing *= np.exp(1j * (-alpha))
#大きさを，中心から右端までが30とする
Wing *= 30 / np.max(np.real(Wing))
#迎角をdegree度にする
Wing *= np.exp(1j * (-degree * rad))


# 直交座標上の位置を得る
#各値を整数に離散化
x_dig = [int(Decimal(i).quantize(Decimal('0'), rounding=ROUND_HALF_UP)) + 201 for i in np.real(Wing)]
y_dig = [int(Decimal(i).quantize(Decimal('0'), rounding=ROUND_HALF_UP)) + 201 for i in np.imag(Wing)]
#重複する点を順序を変えずに除去
xy_list = [(i, j) for i, j in zip(x_dig, y_dig)]
k_pre = ()
for k in xy_list[:]:
    if k == k_pre:
        xy_list.remove(k)
    k_pre = k
x_wall = []
y_wall = []
for xy in xy_list:
    x_wall.append(xy[0])
    y_wall.append(xy[1])
#翼の枠をcsv出力
filename = "./airfoil/wing_data/joukowsky_wall_" + str(degree) + ".csv"
with open(filename, 'w') as file:
    writer = csv.writer(file, lineterminator='\n')
    writer.writerow(x_wall)
    writer.writerow(y_wall)

#翼の枠内の点を追加
multi_check = []
for i in range(len(x_wall)):
    X = x_wall[i]
    if not X in multi_check:
        multi_check.append(X)
        for j in range(len(x_wall)):
            if X == x_wall[j]:
                if y_wall[i] - y_wall[j] > 1:
                    x_wall += [X for k in range(1, y_wall[i] - y_wall[j])]
                    y_wall += [y_wall[j] + k for k in range(1, y_wall[i] - y_wall[j])]
                elif y_wall[j] - y_wall[i] > 1:
                    x_wall += [X for k in range(1, y_wall[j] - y_wall[i])]
                    y_wall += [y_wall[i] + k for k in range(1, y_wall[j] - y_wall[i])]
for i in range(len(y_wall)):
    Y = y_wall[i]
    for j in range(len(y_wall)):
        if Y == y_wall[j]:
            if x_wall[i] - x_wall[j] == 2:
                y_wall += [Y for k in range(1, x_wall[i] - x_wall[j])]
                x_wall += [x_wall[j] + k for k in range(1, x_wall[i] - x_wall[j])]
            elif x_wall[j] - x_wall[i] == 2:
                y_wall += [Y for k in range(1, x_wall[j] - x_wall[i])]
                x_wall += [x_wall[i] + k for k in range(1, x_wall[j] - x_wall[i])]
#重複する点を除去
X = []
Y = []
for xy in list(set([(i, j) for i, j in zip(x_wall, y_wall)])):
    X.append(xy[0])
    Y.append(xy[1])
# 翼形状をcsv出力
filename = "./airfoil/wing_data/joukowsky_" + str(degree) + ".csv"
with open(filename, 'w') as file:
    writer = csv.writer(file, lineterminator='\n')
    writer.writerow(X)
    writer.writerow(Y)

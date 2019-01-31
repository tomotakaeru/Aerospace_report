import numpy as np
import matplotlib.pyplot as plt
import math

#パラメータ
Vw=20
Va=10
alpha=math.pi/3
gamma=40


#格子点を用意
x = np.linspace(-1,1,21)
y = np.linspace(-1,1,21)

X=[]
Y=[]

for i in x:
    for j in y:
        X.append(i)
        Y.append(j)


#速度ベクトルを求める
U=[]
V=[]

for (i,j) in zip(X,Y):
    if i==0 and j==0:
        U.append(0)
        V.append(0)
    else:
        U.append(Vw-Va*math.cos(alpha)-gamma*j/(2*math.pi*(i*i+j*j)))
        V.append(-Va*math.sin(alpha)+gamma*i/(2*math.pi*(i*i+j*j)))

#print(gamma*0.25/(2*math.pi*(0.25*0.25)))


#プロット
plt.figure()
plt.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=90)
plt.plot(0,0,"rx")

# グラフ形式
plt.axes().set_aspect('equal', 'datalim')
plt.xlim([-1.2,1.2])
plt.ylim([-1.2,1.2])
#plt.grid()
plt.draw()
plt.show()
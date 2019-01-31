# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

#初期条件
a=0.6
gamma=5.0
v0=1.0
rad=np.pi/180.0
alpha=5.0*rad
beta=20.0*rad
c=0.5

#計算範囲及び精度
n=181
m=181
dr=0.015
dtheta=2.0*rad
x=np.zeros((n,n))
y=np.zeros((n,n))
X=np.zeros((n,n))
Y=np.zeros((n,n))
Xze=np.zeros((n,n))
Yze=np.zeros((n,n))
stream=np.zeros((n,n))
poten=np.zeros((n,n))

gamma=4.0*np.pi*v0*a*np.sin(alpha+beta) #クッタ条件
Zc=a*np.exp(1j*(np.pi-beta))+c
cgamma=1j*gamma/(2.0*np.pi)


def potential(s,t):
    z=s+1j*t
    f=v0*(z+a**2/z)+cgamma*np.log(z)
    return f

def circle_fc(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alpha)+Zc
    return Z

def Zeta(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alpha)+Zc
    Zeta=Z+c**2/Z
    return Zeta

#ｚ座標系からζ座標系の圧力分布に変換する関数
def pres(s,t):
    z=s+1j*t
    Z=z*np.exp(1j*alpha)+Zc
    f=np.exp(-1j*alpha)*v0*(1-a**2/z**2+1j*gamma/(v0*2*np.pi*z))/(1-c**2/Z**2)
    Cp=1-np.abs(f)**2/v0**2
    return Cp



for i in range(n):
    for k in range(m):
        #z平面上で同心円状に計算点を配置
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=circle_fc(x[i,k],y[i,k])
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

#翼表面部分にi=0が対応
Wing=Zeta(x[0,:],y[0,:])

Cp=pres(x,y)



#翼に働く空気力の計算
cxp=0
cyp=0
Z=circle_fc(x,y)
X=np.real(Z)
Y=np.imag(Z)

Af=0.0
Bf=0.0

#翼表面の圧力分布を積分
for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    Af=Af+Xze[0,l]*fy-Yze[0,l]*fx
    Bf=Bf+fy
    cxp=cxp-cpm*dnx
    cyp=cyp-cpm*dny
    
#抵抗係数、揚力係数
cxp=cxp/(4.0*c)
cyp=cyp/(4.0*c)
cdp=cxp*np.cos(alpha)+cyp*np.sin(alpha)
clp=cyp*np.cos(alpha)-cxp*np.sin(alpha)

print("Cd=",end="")
print(cdp)
print("Cl=",end="")
print(clp)

#風圧中心
xcp=Af/Bf
print("xcp=",end="")
print(xcp)



"""
#翼形状
plt.plot(np.real(Wing),np.imag(Wing),color='k')

#流線
#plt.contour(Xze,Yze,stream,80)

#等圧線
#plt.contour(Xze,Yze,Cp,40)

#圧力分布
plt.plot(Xze[0],Cp[0])

plt.axes().set_aspect('equal', 'datalim')
plt.xlim(-1,1)
plt.ylim(-4,1.5)
#plt.colorbar()
plt.autumn()
plt.show()
"""
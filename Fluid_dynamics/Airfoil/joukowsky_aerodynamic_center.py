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

gamma=4.0*np.pi*v0*a*np.sin(alpha+beta) #Kutta条件
Zc=a*np.exp(1j*(np.pi-beta))+c
cgamma=1j*gamma/(2.0*np.pi)


def potential(s,t):
    z=s+1j*t
    f=v0*(z+a**2/z)+(1j*gamma/(2*np.pi))*np.log(z)
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
        Z=np.exp(1j*alpha)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

Cp=pres(x,y)

Cmc5=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc5[h]=Cmc5[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c**2)



alpha=0*rad
gamma=4.0*np.pi*v0*a*np.sin(alpha+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alpha)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

Cp=pres(x,y)

Cmc0=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc0[h]=Cmc0[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c**2)



alpha=3*rad
gamma=4.0*np.pi*v0*a*np.sin(alpha+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alpha)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

Cp=pres(x,y)

Cmc3=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc3[h]=Cmc3[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c**2)

        

alpha=-3*rad
gamma=4.0*np.pi*v0*a*np.sin(alpha+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alpha)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

Cp=pres(x,y)

Cmc_3=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc_3[h]=Cmc_3[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c**2)

        

alpha=-5*rad
gamma=4.0*np.pi*v0*a*np.sin(alpha+beta)

for i in range(n):
    for k in range(m):
        theta0=0.0
        theta0=-beta
        rrr=a+dr*i
        theta=theta0+dtheta*k
        x[i,k]=rrr*np.cos(theta)
        y[i,k]=rrr*np.sin(theta)
        poten[i,k]=np.real(potential(x[i,k],y[i,k]))
        stream[i,k]=np.imag(potential(x[i,k],y[i,k]))
        Z=np.exp(1j*alpha)*(x[i,k]+1j*y[i,k])+Zc
        Xze[i,k]=np.real(Z+c**2.0/Z)
        Yze[i,k]=np.imag(Z+c**2.0/Z)

Cp=pres(x,y)

Cmc_5=np.zeros(n)
x0=np.linspace(-1.00,0,181)

for l in range(n-1):
    dxw=Xze[0,l+1]-Xze[0,l]
    dyw=Yze[0,l+1]-Yze[0,l]
    dnx=dyw
    dny=-dxw
    cpm=(Cp[0,l+1]+Cp[0,l])/2.0
    fx=-cpm*dnx
    fy=-cpm*dny
    for h in range(n):
        Cmc_5[h]=Cmc_5[h]+((x0[h]-Xze[0,l])*fy+Yze[0,h]*fx)/(0.5*1*1.0**2*16*c**2)


    
plt.plot(x0,Cmc5)
plt.plot(x0,Cmc3)
plt.plot(x0,Cmc0)
plt.plot(x0,Cmc_3)
plt.plot(x0,Cmc_5)
plt.xlim([-0.55,-0.4])
plt.ylim([-1.4,-0.6])
plt.show()

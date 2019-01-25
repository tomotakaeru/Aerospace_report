# -*- coding: utf-8 -*-

from numpy import*
from pylab import*

# 「x」に自分が選んだ機体の値を代入
# 面積、長さ、角度などは自分で計測する必要あり
# 縮尺、rad/degなどに注意

x=1

# 基本諸元＆主翼
h0=3000
U=80.0 #プロペラ機
p=0.909
g=9.8
a0=5.73
CD0=0.04
e=0.8
l=13.36
m=6087
S=28.8
b=17.65
λ=0.412844
Λ=0.03427
Γ=0.1419
c0=2.18
c=2*c0*(1+λ+λ**2)/(3*(1+λ))
AR=b**2/S
aw=a0*cos(Λ)/(1+a0*cos(Λ)/(pi*AR))
CL=2*m*g/(p*U**2*S)
CD=CD0+CL**2/(pi*e*AR)
T=p*U**2*S*CD/2
k1=8/14
k2=13.2/14
caw=0.06
τa=0.3

# 水平尾翼
St=4.52
bt=5.61
λt=0.52
Λt=0.31085
ct=1.2664
c0t=1.583
lt=7.816
Vha=lt*St/(c*S)
ARt=bt**2/St
at=a0*cos(Λt)/(1+a0*cos(Λt)/(pi*ARt))
τe=0.3
Se=1.79

# 垂直尾翼
Sf=4.69
bf=2.28
λf=0.667
ARf=1.55*bf**2/Sf
Λf=0.635
lf=lt
Vfa=lf*Sf/(b*S)
af=a0*cos(Λf)/(1+a0*cos(Λf)/(pi*ARf))
zf=bf*(1+2*λf)/(3*(1+λf))
τr=0.3
Sr=1.40

# その他
εα=2*aw/(pi*AR)
Ixx=0.02*m*b**2
Iyy=0.05*m*l**2
Izz=Ixx+Iyy
Ixz=0.03*Izz
Vfus=23.26
Vfusa=Vfus/(c*S)
h=0.30#
hnw=0.25

# 安定微係数
Cxu=-2*T/(p*S*U**2)-2*CD #プロペラ機
Czα=-aw*(1+at*St*(1-εα)/(aw*S))
Cxα=CL*(1+2*Czα/(pi*e*AR))
Cmα=aw*((h-hnw)-Vha*at*(1-εα)/aw+2*Vfusa/aw)
Cmαdot=-2*Vha*lt*at*εα/c
Czq=-2*Vha*at
Cmq=-2*Vha*lt*at/c
Czδe=-St*at*τe/S
Cmδe=-Vha*at*τe

Cyβ=-Sf*af/S
Clβ=-(1+2*λ)*(aw*Γ+CL*tan(Λ))/(6*(1+λ))
Cnβ=Vfa*af-2*Vfusa*c/b
Clp=-aw*(1+3*λ)/(12*(1+λ))
Cnp=-(1+3*λ)/(12*(1+λ))*Cxα
Cyr=Sf*af*2*lf/(S*b)
Clr=CL/6*(1+3*λ)/(1+λ)+zf*Sf*af*2*lf/(S*b**2)
Cnr=-CD/6*(1+3*λ)/(1+λ)-2*Vfa*af*lf/b
Clδa=aw/6*τa/(1+λ)*(3*(k2**2-k1**2)-2*(1-λ)*(k2**3-k1**3))
Cyδr=Sf*af*τr/S
Clδr=zf/b*Cyδr
Cnδr=-Vfa*af*τr

Xu=p*U*S/(2*m)*Cxu
Zu=p*U*S/(2*m)*(-2*CL)
Xα=p*U**2*S/(2*m)*Cxα
Zα=p*U**2*S/(2*m)*Czα
Mα=p*U**2*S*c/(2*Iyy)*Cmα
Mαdot=p*U*S*c**2/(4*Iyy)*Cmαdot
Zq=p*U*S*c/(4*m)*Czq
Mq=p*U*S*c**2/(4*Iyy)*Cmq
Zδe=p*U**2*S/(2*m)*Czδe
Mδe=p*U**2*S*c/(2*Iyy)*Cmδe

Yβ=p*U**2*S/(2*m)*Cyβ
Lβ=p*U**2*S*b/(2*Ixx)*Clβ
Nβ=p*U**2*S*b/(2*Izz)*Cnβ
Lp=p*U*S*b**2/(4*Ixx)*Clp
Np=p*U*S*b**2/(4*Izz)*Cnp
Yr=p*U*S*b/(4*m)*Cyr
Lr=p*U*S*b**2/(4*Ixx)*Clr
Nr=p*U*S*b**2/(4*Izz)*Cnr
Lδa=p*U**2*S*b/(2*Ixx)*Clδa
Yδr=p*U**2*S/(2*m)*Cyδr
Lδr=p*U**2*S*b/(2*Ixx)*Clδr
Nδr=p*U**2*S*b/(2*Izz)*Cnδr

Lβdash=(Lβ+(Ixz/Ixx)*Nβ)/(1-Ixz**2/(Ixx*Izz))
Nβdash=(Nβ+(Ixz/Izz)*Lβ)/(1-Ixz**2/(Ixx*Izz))
Lpdash=(Lp+(Ixz/Ixx)*Np)/(1-Ixz**2/(Ixx*Izz))
Npdash=(Np+(Ixz/Izz)*Lp)/(1-Ixz**2/(Ixx*Izz))
Lrdash=(Lr+(Ixz/Ixx)*Nr)/(1-Ixz**2/(Ixx*Izz))
Nrdash=(Nr+(Ixz/Izz)*Lr)/(1-Ixz**2/(Ixx*Izz))
Lδadash=Lδa/(1-Ixz**2/(Ixx*Izz))
Nδadash=(Ixz/Izz)*Lδa/(1-Ixz**2/(Ixx*Izz))
Lδrdash=(Lδr+(Ixz/Ixx)*Nδr)/(1-Ixz**2/(Ixx*Izz))
Nδrdash=(Nδr+(Ixz/Izz)*Lδr)/(1-Ixz**2/(Ixx*Izz))



print(CL)
print(CD)
print(T)
print(AR)
print(c)
print(aw)
print("P.2終わり")
print(ARt)
print(Vha)
print(at)
print(ARf)
print(Vfa)
print(af)
print(zf)
print(εα)
print(Ixx)
print(Iyy)
print(Izz)
print(Ixz)
print(Vfus)
print(Vfusa)
print("P.3終わり")
print(Cxu)
print(Cxα)
print(Czα)
print(Cmα)
print(Cmαdot)
print(Czq)
print(Cmq)
print(Czδe)
print(Cmδe)
print("P.4_1/5終わり")
print(Cyβ)
print(Clβ)
print(Cnβ)
print(Clp)
print(Cnp)
print(Cyr)
print(Clr)
print(Cnr)
print(Clδa)
print(Cyδr)
print(Clδr)
print(Cnδr)
print("P.4_2/5終わり")
print(Xu)
print(Zu)
print(Xα)
print(Zα)
print(Mα)
print(Mαdot)
print(Zq)
print(Mq)
print(Zδe)
print(Mδe)
print("P.4_3/5終わり")
print(Yβ)
print(Lβ)
print(Nβ)
print(Lp)
print(Np)
print(Yr)
print(Lr)
print(Nr)
print(Lδa)
print(Yδr)
print(Lδr)
print(Nδr)
print("P.4_4/5終わり")
print(Lβdash)
print(Nβdash)
print(Lpdash)
print(Npdash)
print(Lrdash)
print(Nrdash)
print(Lδadash)
print(Nδadash)
print(Lδrdash)
print(Nδrdash)
print("P.4_5/5終わり")

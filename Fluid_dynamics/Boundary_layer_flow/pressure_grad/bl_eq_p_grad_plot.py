import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(7,4),dpi=100)
ax1 = fig.add_subplot(111)


#諸々厚さ
filepath1='./p_grad/thickness.csv'
df1=pd.read_csv(filepath1)

x1 = [i for i in np.linspace(0,0.5,1001)]
delta = df1.iloc[[0],::50].values[0]
ax1.plot(x1,delta,label="Boundary Layer Thickness",color="b")
delta_star = df1.iloc[[1],::50].values[0]
ax1.plot(x1,delta_star,label="Displacement Thickness",color="m")
theta = df1.iloc[[2],::50].values[0]
ax1.plot(x1,theta,label="Momentum Thickness",color="r")
Theta = df1.iloc[[3],::50].values[0]
ax1.plot(x1,Theta,label="Energy Thickness",color="g")
Cf = df1.iloc[[4],::50].values[0]/2
ax1.plot(x1,Cf,label="Cf",color="k")


#uvベクトル図
X,Y,U,V = [],[],[],[]
x2 = [str(i) for i in np.linspace(0,0.5,11)]
x2[0] = "0"

for i in x2:
    filepath2='./p_grad/uv_x'+ i +'.csv'
    df2=pd.read_csv(filepath2)

    X.extend([i for j in range(len(df2['y'][::20]))])
    Y.extend(df2['y'][::20])
    U.extend( (np.array(df2['u'][::20])/100).tolist() )
    V.extend(df2['v'][::20]/3)

# ax1.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=0.25,color="darkorange")

#網掛け部
yd=Y[0:int(len(X)/len(x2))]
for j in range(len(x2)):
    xd = float(X[j*int(len(yd))])
    xu = [float(X[i])+float(U[i])/5 for i in range( j*int(len(yd)), (j+1)*int(len(yd)) )]
    ax1.fill_betweenx(yd,xu,xd,facecolor="salmon",alpha=0.3)

ax1.set_xlim(0,0.5)
ax1.set_ylim(0,6)
ax1.set_xlabel("$x\;\;$[m]")
ax1.set_ylabel("$y\;\;$[mm]")
ax1.legend(loc="upper right").get_frame().set_alpha(1)

# plt.show()
plt.savefig("hoge2.png",transparent=True)
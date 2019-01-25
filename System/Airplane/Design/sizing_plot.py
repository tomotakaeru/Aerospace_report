import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(7.5,6),dpi=100)
ax1 = fig.add_subplot(111)

x = np.array([i for i in range(0,201)])


CLmaxL = [1.8, 2.4, 3.0] #1.8-3.0
CLmaxTO = [1.6, 2.0, 2.4] #1.6-2.2
color = ["m", "c", "y"]
line = ["-", "--", ":"]

for i in range(len(CLmaxL)): 
    # 着陸性能
    ax1.axvline(1/2*0.125*(101*0.5144)**2*CLmaxL[i]/0.85 *2.2046*0.3048**2, linestyle=line[i], color="k")
for i in range(len(CLmaxTO)): 
    # 離陸性能
    y_to = x*40.3/9300/CLmaxTO[i]
    ax1.plot(x, y_to, color=color[i])
    # 上昇性能
    ax1.axhline(2*((0.033+0.0424*(CLmaxTO[i]/1.2**2)**2)*1.2**2/CLmaxTO[i]+0.024)/0.8, color=color[i], zorder=1)
# 巡航速度
y_cr = [4.78*(4.81/x[i]+x[i]/6366) for i in x]
ax1.plot(x, y_cr, color="g", linestyle="-.")


ax1.quiver(73.2,0.63,-11,0,angles='xy',scale_units='xy',width=0.003,scale=1,color="k")
ax1.quiver(97.6,0.63,-11,0,angles='xy',scale_units='xy',width=0.003,scale=1,color="k")
ax1.quiver(122,0.63,-11,0,angles='xy',scale_units='xy',width=0.003,scale=1,color="k")
ax1.quiver(190,0.515,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="m")
ax1.quiver(190,0.41,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="c")
ax1.quiver(190,0.34,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="y")
ax1.quiver(42,0.58,13,0,angles='xy',scale_units='xy',width=0.003,scale=1,color="g")
ax1.quiver(25,0.25,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="m", zorder=2)
ax1.quiver(20,0.265,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="c", zorder=2)
ax1.quiver(15,0.285,0,0.06,angles='xy',scale_units='xy',width=0.003,scale=1,color="y", zorder=2)

ax1.text(48,0.65,"CLmax,L")
ax1.text(75,0.65,"1.8")
ax1.text(100,0.65,"2.4")
ax1.text(124,0.65,"3.0")
ax1.text(160,0.55,"CLmax,TO")
ax1.text(170,0.485,"1.6")
ax1.text(170,0.39,"2.0")
ax1.text(170,0.325,"2.4")

ax1.plot(121.8, 0.279, "rx" ,markersize="9", markeredgecolor="r", markeredgewidth=2.5)
ax1.text(124.3, 0.295, "${\\bf D}$", fontsize=12)

ax1.set_xlim(0,200)
ax1.set_ylim(0,0.7)
ax1.set_xlabel("$(W/S)_{TO}\;\;[lb/ft^2]$")
ax1.set_ylabel("$(T/W)_{TO}$")

plt.show()
# plt.savefig("sizing_plot.png", transparent=True)
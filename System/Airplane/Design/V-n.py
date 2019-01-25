import numpy as np
import matplotlib.pyplot as plt


# 失速線
rho=1.225
W=1500*9.8
S=16
Va=156/1.944
Vb=122/1.944
x1 = np.arange(0, 9.1, 0.1)
x2 = np.arange(0, 6.1, 0.1)
y1= [1.4*rho*(i*(Va/9))**2/2/W*S for i in x1]
y2= [-1.14*rho*(i*(Vb/6))**2/2/W*S for i in x2]
plt.plot(x1,y1, color='black')
plt.plot(x2,y2, color='black')

# x,y軸
# 矢印の始点
X = 0,-1
Y = -4,0
# 矢印の成分
U = 0,18
V = 11,0
# 矢印
plt.quiver(X,Y,U,V,width=0.004,angles='xy',scale_units='xy',scale=1)

# 直線([x1,x2], [y1,y2])
plt.plot([9,15], [6,6], color='black')
plt.plot([15,15], [6,0], color='black')
plt.plot([6,8], [-3,-3], color='black')
plt.plot([8,15], [-3,0], color='black')

plt.plot([0,9], [6,6], color='black', linestyle='dashed')
plt.plot([9,9], [0,6], color='black', linestyle='dashed')
plt.plot([0,6], [-3,-3], color='black', linestyle='dashed')
plt.plot([6,6], [0,-3], color='black', linestyle='dashed')
plt.plot([8,8], [0,-3], color='black', linestyle='dashed')

# Text挿入
plt.annotate(' 6', xy=(-1,5.8))
plt.annotate('O', xy=(-0.7,-0.5))
plt.annotate('-3', xy=(-1,-3.2))
plt.annotate('A', xy=(8.8,6.3))
plt.annotate('B', xy=(15,6.3))
plt.annotate('C', xy=(7.8,-3.5))
plt.annotate('D', xy=(15,-0.7))
plt.annotate('E', xy=(5.8,-3.5))
plt.annotate('V [ktEAS]', xy=(17.5,-0.2))
plt.annotate('n', xy=(-0.2,7.2))
plt.annotate('LOAD FACTOR', xy=(-1.2,8))
plt.annotate('122', xy=(5.5,0.2))
plt.annotate('152', xy=(7.5,0.2))
plt.annotate('156', xy=(9.3,0.2))
plt.annotate('243', xy=(13.7,0.2))

# 枠線・目盛を消す
plt.axis("off")

# グラフ表示
plt.xlim([-2,19])
plt.ylim([-5,9])
plt.show()

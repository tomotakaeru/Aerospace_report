import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12, 4.5), dpi=100)
fig.subplots_adjust(wspace=0.4)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
N = 50 # 計算する範囲[s]


# 数値解
dt = 0.01
dt_ = 1.0
t1 = np.arange(0, N + 1, dt)
t1_ = np.arange(0, N + 1, dt_)
x1_v = [np.array([[0], [1]])]
x1_v_ = [np.array([[0], [1]])]
x1 = [1]
x1_ = [1]
A = np.array([[0.997, -0.01], [0.01, 1]])
A_ = np.array([[0.364, -0.727], [0.727, 0.582]])
b = np.array([[0.01], [0.]])
b_ = np.array([[0.727], [0.418]])

for i in range(len(t1)-1):
    x1_v.append(np.dot(A, x1_v[i]) + b * np.sin(2 * t1[i]))
    x1.append(float(x1_v[i+1][1]))
for i in range(len(t1_)-1):
    x1_v_.append(np.dot(A_, x1_v_[i]) + b_ * np.sin(2 * t1_[i]))
    x1_.append(float(x1_v_[i+1][1]))

ax1.plot(t1_, x1_, color="m", linestyle="-", label="Numerical Solution (Δt = 1[s])")
ax2.plot(t1, x1, color="m", linestyle="-", label="Numerical Solution (Δt = 0.01[s])")


# 解析解
alpha = -3/20
beta = np.sqrt(391) / 20

t2 = np.arange(0, N + 1, 0.01)
x2 = [-5 / 78 * (np.cos(2 * i) + 5 * np.sin(2 * i)) + 83 / 78 * np.exp(alpha * i) * (np.cos(beta * i) + (alpha + 749 / 830) / beta * np.sin(beta * i)) for i in t2]
ax1.plot(t2, x2, color="c", linestyle="-.", label="Exact Solution")
ax2.plot(t2, x2, color="c", linestyle="-.", label="Exact Solution")



ax1.legend()
ax2.legend()
ax1.set_xlabel("$t$ [s]")
ax2.set_xlabel("$t$ [s]")
ax1.set_ylabel("$x$ [m]")
ax2.set_ylabel("$x$ [m]")
plt.show()
# plt.savefig("MCK_system.png", transparent=True)

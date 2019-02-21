import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipj
from mpl_toolkits.mplot3d import Axes3D


#omega,qの時間遷移
fig11 = plt.figure(figsize=(6, 4), dpi=100)
ax11 = fig11.add_subplot(111)
fig12 = plt.figure(figsize=(6, 4), dpi=100)
ax12 = fig12.add_subplot(111)

filepath1 = "./1/output_omega_q.csv"
df1 = pd.read_csv(filepath1)

omega_x = df1.iloc[:, 0]
omega_y = df1.iloc[:, 1]
omega_z = df1.iloc[:, 2]
q_0 = df1.iloc[:, 3]
q_1 = df1.iloc[:, 4]
q_2 = df1.iloc[:, 5]
q_3 = df1.iloc[:, 6]

dt = 0.01
end_time = 100
t = np.arange(0, end_time + dt, dt)

ax11.plot(t, omega_x, label="$\omega_x$", color="g")
ax11.plot(t, omega_y, label="$\omega_y$", color="m")
ax11.plot(t, omega_z, label="$\omega_z$", color="c")
ax12.plot(t, q_0, label="$q_0$", color="g")
ax12.plot(t, q_1, label="$q_1$", color="m")
ax12.plot(t, q_2, label="$q_2$", color="y")
ax12.plot(t, q_3, label="$q_3$", color="c")

ax11.set_xlim(0, end_time)
ax11.set_ylim(-0.5, 2.0)
ax11.set_xlabel("$Time\;\;$[s]")
ax11.set_ylabel("$\omega\;\;$[rad/s]")
ax11.legend().get_frame().set_alpha(1)

ax12.set_xlim(0, end_time)
ax12.set_ylim(-1, 1)
ax12.set_xlabel("$Time\;\;$[s]")
ax12.set_ylabel("$Quaternion$")
ax12.legend(loc="upper right").get_frame().set_alpha(1)

plt.show()



#omegaの解析解，数値解の比較
fig13 = plt.figure(figsize=(8, 10), dpi=100)
fig13.subplots_adjust(hspace=0.4)
ax131 = fig13.add_subplot(311)
ax132 = fig13.add_subplot(312)
ax133 = fig13.add_subplot(313)

Ix = 1.9
Iy = 1.6
Iz = 2.0
x = 0.1
y = 17.0 * 2.0 * np.pi / 60.0 + 0.1
z = 0.0
E = 0.5 * (Ix * x**2 + Iy * y**2 + Iz * z**2)
L = (Ix * x)**2 + (Iy * y)**2 + (Iz * z)**2
l = np.sqrt((Ix - Iy) * (2 * Iz * E - L) / Ix / Iy / Iz)
k = np.sqrt((Iz - Ix) / (Iy - Ix) * (2 * Iy * E - L) / (2 * Iz * E - L))
a = np.sqrt((2 * Iy * E - L) / Ix / (Iy - Ix))
b = np.sqrt((2 * Iz * E - L) / Iy / (Iz - Iy))
c = np.sqrt((2 * Iy * E - L) / Iz / (Iy - Iz))
K = ellipk(k)
m = k**2

ax131.plot(t, omega_x, color="c", linestyle="-", label="Numerical Solution")
ax131.plot(t, a * ellipj(l * t, m)[1], color="m", linestyle="-.", label="Exact Solution")
ax131.set_xlim(0, end_time)
ax131.set_ylim(-0.1, 0.1)
ax131.set_ylabel("$\omega_x\;\;$[rad/s]")
ax131.legend(loc="upper right").get_frame().set_alpha(1)
ax132.plot(t, omega_y, color="c", linestyle="-", label="Numerical Solution")
ax132.plot(t, b * ellipj(l * t - K, m)[2], color="m", linestyle="-.", label="Exact Solution")
ax132.set_xlim(0, end_time)
ax132.set_ylim(1.8800, 1.8812)
ax132.set_ylabel("$\omega_y\;\;$[rad/s]")
ax132.legend(loc="upper right").get_frame().set_alpha(1)
ax133.plot(t, omega_z, color="c", linestyle="-", label="Numerical Solution")
ax133.plot(t, c * ellipj(l * t, m)[0], color="m", linestyle="-.", label="Exact Solution")
ax133.set_xlim(0, end_time)
ax133.set_ylim(-0.1, 0.1)
ax133.set_xlabel("$Time\;\;$[s]")
ax133.set_ylabel("$\omega_z\;\;$[rad/s]")
ax133.legend(loc="upper right").get_frame().set_alpha(1)

plt.show()



#b-frameの単位ベクトルの軌跡
fig2 = plt.figure(figsize=(11, 8), dpi=100)
fig2.subplots_adjust(wspace=0.1, hspace=0.4)
ax21 = fig2.add_subplot(221, projection="3d")
ax22 = fig2.add_subplot(222, projection="3d")
ax23 = fig2.add_subplot(223, projection="3d")

filepath2 = "./1/output_DCM.csv"
df2 = pd.read_csv(filepath2)

x_x = df2.iloc[:1000:10, 0]
x_y = df2.iloc[:1000:10, 1]
x_z = df2.iloc[:1000:10, 2]
y_x = df2.iloc[:1000:20, 3]
y_y = df2.iloc[:1000:20, 4]
y_z = df2.iloc[:1000:20, 5]
z_x = df2.iloc[:1000:10, 6]
z_y = df2.iloc[:1000:10, 7]
z_z = df2.iloc[:1000:10, 8]

for i in range(len(x_x)):
    ax21.plot([0, x_x[i * 10]], [0, x_y[i * 10]], [0, x_z[i * 10]], "c", linewidth = 0.5)
    ax23.plot([0, z_x[i * 10]], [0, z_y[i * 10]], [0, z_z[i * 10]], "c", linewidth = 0.5)
for i in range(len(y_x)):
    ax22.plot([0, y_x[i * 20]], [0, y_y[i * 20]], [0, y_z[i * 20]], "c", linewidth = 0.5)

ax21.set_xlim(-1, 1)
ax21.set_ylim(-1, 1)
ax21.set_zlim(-1, 1)
ax21.set_xlabel("x")
ax21.set_ylabel("y")
ax21.set_zlabel("z")
ax21.title.set_text("x-axis")
ax21.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax21.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax21.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax21.view_init(elev=25, azim=30)
ax22.set_xlim(-1, 1)
ax22.set_ylim(-1, 1)
ax22.set_zlim(-1, 1)
ax22.set_xlabel("x")
ax22.set_ylabel("y")
ax22.set_zlabel("z")
ax22.title.set_text("y-axis")
ax22.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax22.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax22.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax22.view_init(elev=25, azim=30)
ax23.set_xlim(-1, 1)
ax23.set_ylim(-1, 1)
ax23.set_zlim(-1, 1)
ax23.set_xlabel("x")
ax23.set_ylabel("y")
ax23.set_zlabel("z")
ax23.title.set_text("z-axis")
ax23.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax23.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax23.set_zticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax23.view_init(elev=25, azim=30)

plt.show()

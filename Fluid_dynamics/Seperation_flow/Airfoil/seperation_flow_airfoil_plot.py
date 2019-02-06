import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
from matplotlib.animation import FuncAnimation, PillowWriter


# 計算格子の用意
x_min = -11.0
x_max = 33.0
y_min = -11.0
y_max = 11.0
n_last = 5000
x_np = np.linspace(x_min, x_max, 81)
y_np = np.linspace(y_max, y_min, 41)
p_np = []
x, y = np.meshgrid(x_np, y_np)
X = []
Y = []
for j in y_np:
    for i in x_np:
        X.append(i)
        Y.append(j)


def read_p_csv(file_name):
    p_lis = []
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            p_row = []
            for p in row:
                p_row.append(round(float(p), 3))
            p_lis.append(p_row)
    p_np_part = np.array(p_lis)
    return p_np_part

def read_u_csv(file_name):
    u_lis = []
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        for index_y, row in enumerate(reader):
            if index_y%10 ==0:
                u_row = []
                for index_x, u in enumerate(row):
                    if index_x%10 == 0:
                        u_row.append(round(float(u), 3))
                u_lis.append(u_row)
    u_np_part = np.array(u_lis)
    return u_np_part

def read_wing_wall_csv(file_name):
    wing_list = []
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            row_list = []
            for i in row:
                row_list.append(int(i))
            wing_list.append(np.array(row_list))
    wing_wall_x = wing_list[0] * (x_max - x_min) / 801 + x_min
    wing_wall_y = wing_list[1] * (y_max - y_min) / 401 + y_min
    return wing_wall_x, wing_wall_y


def plot_p_gif(alpha):
    fig = plt.figure()
    ims = []
    for i in range(1, 200):
        print(i)
        p_np_part = read_p_csv("./airfoil/output_csv/p_" + alpha + "_" + str(i * 50) + ".csv")
        p_np.append(p_np_part)
        im = plt.imshow(p_np[-1], extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.2, vmax=0.15, cmap="Spectral_r")  # vmin,vmaxで等高線のカラーリングを調整
        ims.append([im])
    anim = animation.ArtistAnimation(fig, ims, interval=30, blit=False)
    anim.save("./airfoil/output_image/karman_" + alpha + ".gif", writer=PillowWriter())
    plt.show()

def plot_iso_pressure(alpha, time):
    p = read_p_csv("./airfoil/output_csv/p_" + alpha + "_" + time + ".csv")
    plt.imshow(p, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.2, vmax=0.15, cmap="Spectral_r")  # vmin,vmaxで等高線のカラーリングを調整
    plt.gca().set_aspect('equal')
    plt.show()

def plot_streamline(alpha, time):
    wing_wall_x, wing_wall_y = read_wing_wall_csv("./airfoil/wing_data/joukowsky_wall_" + alpha + ".csv")
    U = read_u_csv("./airfoil/output_csv/u_" + alpha + "_" + time + ".csv")
    V = np.array(read_u_csv("./airfoil/output_csv/v_" + alpha + "_" + time + ".csv")) * 2
    plt.figure()
    plt.plot(wing_wall_x, wing_wall_y)
    plt.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=2, width=0.001)
    plt.xlim(-5, 30)
    plt.ylim(-5, 5)
    plt.axes().set_aspect('equal')
    plt.show()

def plot_cdclcp(alpha):
    cd, cl, cp1, cp2 = [], [], [], []
    t = np.linspace(1, 10000, 10000)
    with open("./airfoil/output_csv/cdclcp_" + alpha + ".csv", 'r') as f:
        reader = csv.reader(f)
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            cd.append(float(row[0]))
            cl.append(float(row[1]))
            cp1.append(float(row[2]))
            # cp2.append(float(row[3]))
    plt.plot(t, cd, label="Cd")
    plt.plot(t, cl, label="Cl")
    plt.plot(t, cp1, label="Cp")
    # plt.plot(t, cp2, label="Cp2")
    plt.xlabel("step")
    plt.xlim([0, 5000])
    plt.ylim([-1, 4])
    plt.legend()
    plt.show()




if __name__ == "__main__":
    alpha = "20"  #適切な迎角
    time = "8000"  #切り取りたい時刻

    plot_p_gif(alpha)
    plot_iso_pressure(alpha, time)
    plot_streamline(alpha, time)
    plot_cdclcp(alpha)

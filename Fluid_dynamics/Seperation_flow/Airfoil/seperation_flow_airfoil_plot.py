import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
from matplotlib.animation import FuncAnimation, PillowWriter


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
        for index_y, row in enumerate(reader):  # y座標を固定したときの各xの圧力pのリストがrow
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
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            row_list = []
            for i in row:
                row_list.append(int(i))
            wing_list.append(np.array(row_list))
    wing_wall_x = wing_list[0] * (x_max - x_min) / 801 + x_min
    wing_wall_y = wing_list[1] * (y_max - y_min) / 401 + y_min
    return wing_wall_x, wing_wall_y


def plot_p_gif():
    fig = plt.figure()
    ims = []
    for i in range(1, 200):
        print(i)
        p_np_part = read_p_csv("./airfoil/data/p_10_" + str(i * 50) + ".csv")
        p_np.append(p_np_part)
        im = plt.imshow(p_np[-1], extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.2, vmax=0.2, cmap="Spectral_r")  # vmin,vmaxで等高線のカラーリングを調整
        ims.append([im])
    anim = animation.ArtistAnimation(fig, ims, interval=30, blit=False)
    # anim.save('/airfoil/karman_wing_10.gif', writer=PillowWriter)
    plt.show()

def plot_iso_pressure():
    p = read_p_csv("./airfoil/data/p_10_7000.csv")
    plt.imshow(p, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.2, vmax=0.2, cmap="Spectral_r")  # vmin,vmaxで等高線のカラーリングを調整
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.show()

def plot_streamline():
    wing_wall_x, wing_wall_y = read_wing_wall_csv("./airfoil/joukowsky_wall_10.csv")
    U = read_u_csv("./airfoil/data/u_10_50.csv")
    V = np.array(read_u_csv("./airfoil/data/v_10_50.csv")) * 2
    plt.figure()
    plt.plot(wing_wall_x, wing_wall_y)
    plt.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=2, width=0.001)
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()

def plot_cdclcp():
    cd, cl, cp1, cp2 = [], [], [], []
    t = np.linspace(1, 10000, 10000)
    with open("./airfoil/data/cdclcp_10.csv", 'r') as f:
        reader = csv.reader(f)
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            cd.append(float(row[0]))
            cl.append(float(row[1]))
            cp1.append(float(row[2]))
            cp2.append(float(row[3]))
    plt.plot(t, cd, label="Cd")
    plt.plot(t, cl, label="Cl")
    plt.plot(t, cp1, label="Cp1")
    plt.plot(t, cp2, label="Cp2")
    plt.xlabel("step")
    plt.xlim([0, 5000])
    plt.ylim([-1, 3])
    plt.legend()
    plt.show()




if __name__ == "__main__":
    plot_p_gif()
    plot_iso_pressure()
    plot_streamline()
    plot_cdclcp()

import numpy as np
import matplotlib.pyplot as plt


#条件設定
rho = 2700  #密度[kg/m^3] 
E = 70 * 1e9  #ヤング率[Pa]
h = 0.002  #梁の厚さ[m]
l = 1  #梁の長さ[m]
lambda_i = [4.7300408, 7.8532046, 10.9956078]  #coshλcosλ=1の解（一様断面はりの流用）
deg = len(lambda_i)  #モードの次数
M = np.zeros((deg, deg))  #質量行列の準備
K = np.zeros((deg, deg))  #剛性行列の準備


#テーパー梁の厚さ
def h_tp(x):
    return h * (1 + 2 * x / l)

#両端自由の境界条件を満たす関数φj
def calc_phi(i, x):
    try:
        lambda_j = lambda_i[i]
        a_j = lambda_j / l
    except:
        print("iは0-2に限定されてるよ")
    else:
        return (np.cosh(a_j * x) + np.cos(a_j * x)) - (np.sinh(lambda_j) + np.sin(lambda_j)) / (np.cosh(lambda_j) - np.cos(lambda_j)) * (np.sinh(a_j * x) + np.sin(a_j * x))

#φjの2階微分
def calc_ddphi(i, x):
    try:
        lambda_j = lambda_i[i]
        a_j = lambda_j / l
    except:
        print("iは0-2に限定されてるよ")
    else:
        return a_j ** 2 * ((np.cosh(a_j * x) - np.cos(a_j * x)) - (np.sinh(lambda_j) + np.sin(lambda_j)) / (np.cosh(lambda_j) - np.cos(lambda_j)) * (np.sinh(a_j * x) - np.sin(a_j * x)))

#質量行列を計算
def calc_M(dx=0.001):
    for i in range(deg):
        for j in range(deg):
            m=0
            for n in range(int(l / dx)):  #各要素を台形近似で積分
                x_mid = (n + 0.5) * dx  #各台形の中央位置
                m += rho * h_tp(x_mid) * calc_phi(i, x_mid) * calc_phi(j, x_mid) * dx
            M[i, j] = m
    return M

#剛性行列を計算
def calc_K(dx=0.001):
    for i in range(deg):
        for j in range(deg):
            k=0
            for n in range(int(l / dx)):  #各要素を台形近似で積分
                x_mid = (n + 0.5) * dx  #各台形の中央位置
                k += E / 12 * h_tp(x_mid)** 3 * calc_ddphi(i, x_mid) * calc_ddphi(j, x_mid) * dx
            K[i, j] = k
    return K

#テーパー梁の固有振動数ωの2乗とモードを出力
def calc_eigen_tp():
    return np.linalg.eig(np.dot(np.linalg.inv(calc_M()), calc_K()))

#固有振動モードを求める
def mode1_tp(x): #1次モード
    return sum([C[0, i] * calc_phi(i, x) for i in range(deg)])
def mode2_tp(x): #2次モード
    return sum([C[1, i] * calc_phi(i, x) for i in range(deg)])
def mode3_tp(x): #3次モード
    return sum([C[2, i] * calc_phi(i, x) for i in range(deg)])

#テーパー梁の1-3次モードをプロット
def plot_tp():
    fig = plt.figure(figsize=(8, 7), dpi=100)
    ax1 = fig.add_subplot(111)

    X = np.arange(0, l, 0.001)
    ax1.plot(X,-mode1_tp(X), label="mode1: $\omega1=$" + str(int(omega_tp[0])) + " [rad/s]")
    ax1.plot(X, mode2_tp(X), label="mode2: $\omega2=$" + str(int(omega_tp[1])) + " [rad/s]")
    ax1.plot(X, mode3_tp(X), label="mode3: $\omega3=$" + str(int(omega_tp[2])) + " [rad/s]")
    ax1.set_title("Tapered Beam Mode")
    ax1.set_xlabel("$x\;$[m]")
    ax1.set_ylabel("$z/C\;$[m]")
    ax1.legend(loc="upper center").get_frame().set_alpha(1)
    plt.show()

#一様断面梁の1-3次モードをプロット
def plot_uni():
    fig = plt.figure(figsize=(8, 7), dpi=100)
    ax1 = fig.add_subplot(111)

    X = np.arange(0, l, 0.001)
    ax1.plot(X, calc_phi(0, X), label="mode1: $\omega1=$" + str(int(omega_uni[0])) + " [rad/s]")
    ax1.plot(X, calc_phi(1, X), label="mode2: $\omega2=$" + str(int(omega_uni[1])) + " [rad/s]")
    ax1.plot(X, calc_phi(2, X), label="mode3: $\omega3=$" + str(int(omega_uni[2])) + " [rad/s]")
    ax1.set_title("Uniform Beam Mode")
    ax1.set_xlabel("$x\;$[m]")
    ax1.set_ylabel("$z/C\;$[m]")
    ax1.legend(loc="upper center").get_frame().set_alpha(1)
    plt.show()

#モードの比較をプロット
def plot_compare_tp_uni():
    fig = plt.figure(figsize=(9, 9), dpi=100)
    fig.subplots_adjust(wspace=0.3, hspace=0.5)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)

    X = np.arange(0, l, 0.001)

    ax1.plot(X, -mode1_tp(X), label="Tapered: " + str(int(omega_tp[0])) + " [rad/s]")
    ax1.plot(X, calc_phi(0, X), label="Uniform: " + str(int(omega_uni[0])) + " [rad/s]")
    ax1.set_title("Mode 1")
    ax1.set_xlabel("$x\;$[m]")
    ax1.set_ylabel("$z/C\;$[m]")
    ax1.legend(loc="upper center", prop={'size':7}).get_frame().set_alpha(1)

    ax2.plot(X, mode2_tp(X), label="Tapered: " + str(int(omega_tp[1])) + " [rad/s]")
    ax2.plot(X, calc_phi(1, X), label="Uniform: " + str(int(omega_uni[1])) + " [rad/s]")
    ax2.set_title("Mode 2")
    ax2.set_xlabel("$x\;$[m]")
    ax2.set_ylabel("$z/C\;$[m]")
    ax2.legend(loc="upper center", prop={'size':7}).get_frame().set_alpha(1)

    ax3.plot(X, mode3_tp(X), label="Tapered: " + str(int(omega_tp[2])) + " [rad/s]")
    ax3.plot(X, calc_phi(2, X), label="Uniform: " + str(int(omega_uni[2])) + " [rad/s]")
    ax3.set_title("Mode 3")
    ax3.set_xlabel("$x\;$[m]")
    ax3.set_ylabel("$z/C\;$[m]")
    ax3.legend(loc="upper center", prop={'size':7}).get_frame().set_alpha(1)

    plt.show()


if __name__ == "__main__":
    omega2, C = calc_eigen_tp()
    omega_tp = np.sqrt(omega2)
    omega_uni = [(lambda_i[i] / l)** 2 * np.sqrt(E * h ** 2 / 3 / rho) for i in range(deg)]
    print(omega_tp, omega_uni)

    plot_tp()
    plot_uni()
    plot_compare_tp_uni()

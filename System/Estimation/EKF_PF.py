import numpy as np
import matplotlib.pyplot as plt


#状況設定
rho_0 = 1.225  #[kg/m^3]
g = 9.807  #[m/s^2]
D = 100  #[kg/m^2]
beta = 1 / 7000  #[1/m]
h_0 = 200e3  #[m]
V_reentry = 12e3  #[m/s]
gamma = -5 / 180 * np.pi  #[rad]
Q = np.diag([10, 10, 5, 5])
R = np.diag([100, 1 / 180 * np.pi])
m_0 = np.array([0, V_reentry * np.cos(gamma), h_0, V_reentry * np.sin(gamma)])
V_0 = np.zeros((4, 4))
spot_obs = 1350e3  #[m]

#計算の準備
x_deg = 4  #状態量の次元
y_deg = 2  #観測量の次元
dt = 1  #[s]
end_time = 1000  #[s]
time = np.arange(0, end_time, dt)
N = int(end_time / dt + 1)  #データ数
X_act = np.zeros((N, x_deg))  #実際の状態量の入れ物
Y_act = np.zeros((N, y_deg))  #観測量の入れ物
X_act[0] = m_0

#システム
def F(x_t):
    tmp1 = rho_0 / 2 / D * np.e ** (-beta * x_t[2]) * x_t[1]** 2
    tmp2 = rho_0 / 2 / D * np.e ** (-beta * x_t[2]) * x_t[3]** 2
    return np.array([x_t[0] + x_t[1] * dt, x_t[1] + -np.sign(x_t[1]) * tmp1 * dt, x_t[2] + x_t[3] * dt, x_t[3] + (-np.sign(x_t[3]) * tmp2 - g) * dt])
def G(x_t):
    return np.array([np.sqrt((x_t[0] - spot_obs)** 2 + x_t[2]** 2), np.arctan(x_t[2] / (spot_obs - x_t[0]))])

#実際の状態量を生成
def simulate(X,Y,Q,R,N):
    Y[0] = G(X[0]) + np.array([np.random.normal(0, i) for i in np.diag(R)])
    for i in range(N - 1):
        X[i + 1] = F(X[i]) + np.array([np.random.normal(0, i) for i in np.diag(Q)])
        Y[i + 1] = G(X[i + 1]) + np.array([np.random.normal(0, i) for i in np.diag(R)])

        if X[i + 1][2] < 0:  #地表に着いたら計算を止める
            X[i + 1][1:] = np.zeros(x_deg - 1)
            return np.delete(X, np.s_[(i + 2):N], 0), np.delete(Y, np.s_[(i + 2):N], 0)

#システムの偏微分を計算
def calc_partial_F(m_t):
	return np.array([[1, dt, 0, 0],
	 [0, 1 - rho_0 / D * np.e ** (-beta * m_t[2]) * m_t[1] * dt, rho_0 / 2 / D * beta * np.e ** (-beta * m_t[2]) * m_t[1] ** 2 * dt, 0],
	 [0, 0, 1, dt],
	 [0, rho_0 / D * np.e ** (-beta * m_t[2]) * m_t[1] * dt, -rho_0 / 2 / D * beta * np.e ** (-beta * m_t[2]) * m_t[1] ** 2 * dt, 1]])
def calc_partial_G(m_t):
	return np.array([[-(spot_obs - m_t[0]) / np.sqrt((spot_obs - m_t[0])** 2 + m_t[2]** 2), 0, m_t[2] / np.sqrt((spot_obs - m_t[0])** 2 + m_t[2]** 2), 0],
         [m_t[2] / ((spot_obs - m_t[0])** 2) / (1 + (m_t[2] / (spot_obs - m_t[0]))** 2), 0, 1 / (spot_obs - m_t[0]) / (1 + (m_t[2] / (spot_obs - m_t[0]))** 2), 0]])

#EKF
def estimate_EKF(m, V, Q, R, N, Y_act):
    for i in range(N - 1):
        m[i + 1] = F(m[i])
        A = calc_partial_F(m[i])
        V[i + 1] = np.dot(np.dot(A, V[i]), A.T) + Q

        C = calc_partial_G(m[i + 1])
        K = np.dot(np.dot(V[i + 1], C.T), np.linalg.inv(np.dot(np.dot(C, V[i + 1]), C.T) + R))
        m[i + 1] = m[i + 1] + np.dot(K, Y_act[i].T - G(m[i + 1]))
        V[i + 1] = np.dot((np.eye(x_deg) - np.dot(K, C)), V[i + 1])
    return m, V

#PF
def estimate_PF(X, w, Q, R, N, Y_act):
    for k in range(N - 1):
        for j in range(N_particle):
            X[k + 1][j] = F(X[k][j]) + np.array([np.random.normal(0, i) for i in np.diag(Q)])
            w[k + 1][j] = np.clip(w[k][j] * np.exp(-sum([(Y_act[k + 1][i] - G(X[k + 1][j])[i])** 2 / np.diag(R)[i] for i in range(y_deg)])), 1e-300, 100)
        w_sum = sum(w[k + 1])
        for j in range(N_particle):
            w[k + 1][j] /= w_sum
    return X, w



if __name__ == "__main__":
    #実際の状態量を生成
    X_act, Y_act = simulate(X_act, Y_act, Q, R, N)
    N_act = len(X_act)
    print("到着までの時間: {0:}[s]".format((N_act - 1) * dt))

    x_act = [X_act[i][0] / 1000 for i in range(N_act)]
    h_act = [X_act[i][2] / 1000 for i in range(N_act)]
    x_obs = [spot_obs / 1000 - Y_act[i][0] / 1000 * np.cos(Y_act[i][1]) for i in range(len(Y_act))]
    h_obs = [Y_act[i][0] / 1000 * np.sin(Y_act[i][1]) for i in range(len(Y_act))]


    #EKFの準備
    m_est = np.zeros((N_act, x_deg))  #平均値の入れ物
    V_est = np.zeros((N_act, x_deg, x_deg))  #共分散行列の入れ物
    m_est[0] = m_0
    V_est[0] = V_0

    #EKFによる推定
    m_est, V_est = estimate_EKF(m_est, V_est, Q, R, N_act, Y_act)

    x_est_EKF = [m_est[i][0] / 1000 for i in range(N_act)]
    xdot_est_EKF = [m_est[i][1] / 1000 for i in range(N_act)]
    h_est_EKF = [m_est[i][2] / 1000 for i in range(N_act)]
    hdot_est_EKF = [m_est[i][3] / 1000 for i in range(N_act)]


    #PFの準備
    N_particle = 1000
    X_est = np.zeros((N_act, N_particle, x_deg))
    w = np.zeros((N_act, N_particle))
    X_est[0] = np.array([m_0 for i in range(N_particle)])
    w[0] = np.array([1/N_particle for i in range(N_particle)])

    #PFによる推定
    X_est, w = estimate_PF(X_est, w, Q, R, N_act, Y_act)

    x_est_PF = [m_0[0] / 1000] + [sum([w[k][j] * X_est[k][j][0] for j in range(N_particle)]) / 1000 for k in range(1, N_act)]
    xdot_est_PF = [m_0[1] / 1000] + [sum([w[k][j] * X_est[k][j][1] for j in range(N_particle)]) / 1000 for k in range(1, N_act)]
    h_est_PF = [m_0[2] / 1000] + [sum([w[k][j] * X_est[k][j][2] for j in range(N_particle)]) / 1000 for k in range(1, N_act)]
    hdot_est_PF = [m_0[3] / 1000] + [sum([w[k][j] * X_est[k][j][3] for j in range(N_particle)]) / 1000 for k in range(1, N_act)]


    #EKF,PFの性能比較(平均二乗誤差)
    e_EKF = np.mean([(x_est_EKF[i] - x_act[i])** 2 + (h_est_EKF[i] - h_act[i])** 2 for i in range(N_act)])
    e_PF  = np.mean([(x_est_PF[i]  - x_act[i])** 2 + (h_est_PF[i]  - h_act[i])** 2 for i in range(N_act)])
    print("実際の軌道との平均二乗誤差")
    print("EKFの場合: {0}, PFの場合: {1}".format(e_EKF, e_PF))


    #軌跡をプロット
    fig = plt.figure(figsize=(15, 4), dpi=100)
    ax1 = fig.add_subplot(111)

    ax1.plot(x_act, h_act, label="Actual", color="m", zorder=2)
    ax1.plot(x_obs, h_obs, label="Observed", color="c", zorder=1)
    ax1.plot(x_est_EKF, h_est_EKF, label="EKF", color="y", zorder=4)
    ax1.plot(x_est_PF, h_est_PF, label="PF", color="r", zorder=3)

    ax1.scatter(1350, 0, s=150, marker="+", color="g")  #観測地点
    ax1.scatter(x_act[-1], 0, s=150, marker="x", color="k", zorder=5)  #到着地点
    ax1.set_xlim(-10, 1400)
    ax1.set_ylim(-5, 205)
    ax1.set_xlabel("$x\;[km]$")
    ax1.set_ylabel("$h\;[km]$")
    ax1.legend(loc="upper right").get_frame().set_alpha(1)
    plt.show()


    #時系列データをプロット
    fig = plt.figure(figsize=(11, 8), dpi=100)
    fig.subplots_adjust(wspace=0.5, hspace=0.3)
    Time = np.linspace(0, (N_act - 1) * dt, N_act)

    #t-x
    ax1 = fig.add_subplot(221)
    ax1.plot(Time, x_est_EKF, label="EKF", color="y", zorder=4)
    ax1.plot(Time, x_est_PF, label="PF", color="r", zorder=3)

    ax1.set_xlim(0, 400)
    ax1.set_ylim(0, 1400)
    ax1.set_xlabel("$Time\;[s]$")
    ax1.set_ylabel("$x\;[km]$")
    ax1.legend(loc="upper left").get_frame().set_alpha(1)

    #t-x_dot
    ax2 = fig.add_subplot(222)
    ax2.plot(Time, xdot_est_EKF, label="EKF", color="y", zorder=4)
    ax2.plot(Time, xdot_est_PF, label="PF", color="r", zorder=3)

    ax2.set_xlim(0, 400)
    ax2.set_ylim(0, 13)
    ax2.set_xlabel("$Time\;[s]$")
    ax2.set_ylabel("$x\_dot\;[km/s]$")
    ax2.legend(loc="upper right").get_frame().set_alpha(1)

    #t-h
    ax3 = fig.add_subplot(223)
    ax3.plot(Time, h_est_EKF, label="EKF", color="y", zorder=4)
    ax3.plot(Time, h_est_PF, label="PF", color="r", zorder=3)

    ax3.set_xlim(0, 400)
    ax3.set_ylim(0, 200)
    ax.set_xlabel("$Time\;[s]$")
    ax3.set_ylabel("$h\;[km]$")
    ax3.legend(loc="upper right").get_frame().set_alpha(1)
    
    #t-h_dot
    ax4 = fig.add_subplot(224)
    ax4.plot(Time, hdot_est_EKF, label="EKF", color="y", zorder=4)
    ax4.plot(Time, hdot_est_PF, label="PF", color="r", zorder=3)

    ax4.set_xlim(0, 400)
    ax4.set_ylim(-2, 0.2)
    ax4.set_xlabel("$Time\;[s]$")
    ax4.set_ylabel("$h\_dot\;[km/s]$")
    ax4.legend(loc="upper left").get_frame().set_alpha(1)

    
    plt.show()

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize


#観測データ取得
filepath = './Hubbles_constant.csv'
df = pd.read_csv(filepath)
Distance = df["distance"]
Velocity = df["recession_velocity"]


#最小二乗法でパラメータ推定
def calc_param_no_const(X, Y):
    return [sum(X * Y) / sum(X * X)]

a_no_const = calc_param_no_const(Distance, Velocity)

def calc_param_have_const(X, Y):
    one = np.array([[1 for i in range(len(X))]])
    X = np.array([list(X)])
    Y = np.array(Y).T
    X = np.vstack((one, X)).T
    return np.linalg.solve(np.dot(X.T, X), np.dot(X.T, Y))

a_have_const = calc_param_have_const(Distance, Velocity)


#推定したパラメータでのフィッティングデータ
Distance_model = np.arange(0, 2.1, 0.001)
Velocity_model_no_const = [a_no_const[0] * Distance_model[i] for i in range(len(Distance_model))]
Velocity_model_have_const = [a_have_const[0] + a_have_const[1] * Distance_model[i] for i in range(len(Distance_model))]


#プロット
fig = plt.figure(figsize=(7, 5.5), dpi=100)
ax1 = fig.add_subplot(111)

ax1.scatter(Distance, Velocity, label="Data", color="m", marker="x")
ax1.plot(Distance_model, Velocity_model_no_const, label="Estimation (no constant)", color="c")
ax1.plot(Distance_model, Velocity_model_have_const, label="Estimation (have constant)", color="y")

ax1.set_xlim(0, 2.1)
ax1.set_ylim(-400, 1200)
ax1.set_xlabel("$Distance\;\;$[Mpc]")
ax1.set_ylabel("$Recession\;Velocity\;\;$[km/s]")
ax1.legend(loc="upper left").get_frame().set_alpha(1)
plt.show()


#AIC（赤池情報量基準）によりモデルを評価する
def calc_AIC(X, Y, param):
    n = len(Y)
    if len(param) == 1:
        return n * np.log(sum([(Y[i] - param[0] * X[i])** 2 for i in range(n)]) / n) + 2 * 1
    elif len(param) == 2:
        return n * np.log(sum([(Y[i] - param[0] - param[1] * X[i])** 2 for i in range(n)]) / n) + 2 * 2
    else:
        print("ERROR: can't accept model except *Simple Linear Regression*")
def calc_AICc(X, Y, param):
    n = len(Y)
    if len(param) == 1:
        return n * np.log(sum([(Y[i] - param[0] * X[i])** 2 for i in range(n)]) / n) + 2 * 1 * n / (n - 2 * 1 - 1)
    elif len(param) == 2:
        return n * np.log(sum([(Y[i] - param[0] - param[1] * X[i])** 2 for i in range(n)]) / n) + 2 * 2 * n / (n - 2 * 2 - 1)
    else:
        print("ERROR: can't accept model except *Simple Linear Regression*")

AIC_no_const = calc_AIC(Distance, Velocity, a_no_const)
AIC_have_const = calc_AIC(Distance, Velocity, a_have_const)
AICc_no_const = calc_AICc(Distance, Velocity, a_no_const)
AICc_have_const = calc_AICc(Distance, Velocity, a_have_const)


#print
print(a_no_const, a_have_const)
print(sum([(Velocity[i] - a_no_const[0] * Distance[i])** 2 for i in range(len(Velocity))]), sum([(Velocity[i] - a_have_const[0] - a_have_const[1] * Distance[i])** 2 for i in range(len(Velocity))])) #残差二乗和
print(AIC_no_const, AIC_have_const, AICc_no_const, AICc_have_const)

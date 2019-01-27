import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize


#観測データ取得
filepath = './satellite-temp.csv'
df = pd.read_csv(filepath)
Time = df.iloc[:, 0]
Temp = df.iloc[:, 1]


#各データとパラメータから確率を計算
def calc_probability(param, time_i, data_i):
    e = data_i - param[0] - param[1] * np.sin(param[2] * time_i + param[3])
    return np.exp(-e ** 2 / 2 / param[4] ** 2) / np.sqrt(2 * np.pi * param[4] ** 2)

#全データについて対数尤度を計算
def calc_LL(param, time, data):
    LL = 0
    for i in range(len(time)):
        LL += np.log(calc_probability(param, time[i], data[i]))
    return -LL  #minimizeによって求めるため

#最尤推定
param_init = [12, 3, 2 * np.pi / 6000, 0, 1]
param_bound = [(10, 15), (1, 5), (2 * np.pi / 10000, 2 * np.pi / 3000), (0, 2 * np.pi), (0, None)]
param_MLE = optimize.minimize(calc_LL, param_init, (Time, Temp), method='L-BFGS-B', bounds=param_bound).x
print("T0: {0[0]:.1f}, a: {0[1]:.2f}, ω: {0[2]:.2f}, θ: {0[3]:.2f}, σ: {0[4]:.2f}".format(param_MLE))

#推定したパラメータでのフィッティングデータ
Time_model = np.arange(14000)
Temp_model = [param_MLE[0] + param_MLE[1] * np.sin(param_MLE[2] * Time_model[i] + param_MLE[3]) for i in range(len(Time_model))]


#プロット
fig = plt.figure(figsize=(8, 5), dpi=100)
ax1 = fig.add_subplot(111)

ax1.scatter(Time, Temp, label="Satellite Temperature Data", color="m", marker="x")
ax1.plot(Time_model, Temp_model, label="Satellite Temperature Estimation", color="c")

ax1.set_xlim(0, 14000)
ax1.set_ylim(8, 18)
ax1.set_xlabel("$Time\;\;$[s]")
ax1.set_ylabel("$Temperature\;\;$[deg]")
ax1.legend(loc="upper right").get_frame().set_alpha(1)
plt.show()

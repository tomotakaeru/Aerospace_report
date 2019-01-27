import numpy as np
import matplotlib.pyplot as plt


#データ
x = [-2, -1, 0, 1, 2]
y = [-1.8623, 0.6339, -2.2588, 2.0622, 2.7188]
noise_distribution = [0, 1] #[mean, variance]
apriori_distribution = [1, 0.09]

#推定(LSM,MAP)
def calc_param_LSM(x, y):
    return sum([y[i] * x[i] for i in range(len(x))]) / sum([x[i]** 2 for i in range(len(x))])
a_LSM = calc_param_LSM(x, y)

def calc_param_MAP(x, y, noise, apriori):
    return (sum([y[i] * x[i] for i in range(len(x))]) / noise[1] + apriori[0] / apriori[1]) / (sum([x[i]** 2 for i in range(len(x))]) / noise[1] + 1 / apriori[1])
a_MAP = calc_param_MAP(x, y, noise_distribution, apriori_distribution)

#推定したパラメータでのフィッティングデータ
x_model = np.arange(-2.1, 2.1, 0.01)
y_model_LSM = [a_LSM * x_model[i] for i in range(len(x_model))]
y_model_MAP = [a_MAP * x_model[i] for i in range(len(x_model))]

#プロット
fig = plt.figure(figsize=(6, 5), dpi=100)
ax1 = fig.add_subplot(111)

ax1.scatter(x, y, label="Data", color="m", marker="x")
ax1.plot(x_model, y_model_LSM, label="Estimation (MLE)", color="c")
ax1.plot(x_model, y_model_MAP, label="Estimation (MAP)", color="y")

ax1.set_xlim(-2.1, 2.1)
ax1.set_ylim(-3, 3)
ax1.set_xlabel("$x$")
ax1.set_ylabel("$y$")
ax1.legend(loc="upper left").get_frame().set_alpha(1)
plt.show()

print(a_LSM, a_MAP)

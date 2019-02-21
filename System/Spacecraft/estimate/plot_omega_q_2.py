import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#omega,omega_hat,domega,q,q_hat,dqの時間遷移
fig11 = plt.figure(figsize=(8, 6), dpi=100)
fig11.subplots_adjust(hspace=0.4)
ax111 = fig11.add_subplot(211)
ax112 = fig11.add_subplot(212)
fig12 = plt.figure(figsize=(8, 6), dpi=100)
fig12.subplots_adjust(hspace=0.4)
ax121 = fig12.add_subplot(211)
ax122 = fig12.add_subplot(212)
fig13 = plt.figure(figsize=(8, 6), dpi=100)
fig13.subplots_adjust(hspace=0.4)
ax131 = fig13.add_subplot(211)
ax132 = fig13.add_subplot(212)
fig21 = plt.figure(figsize=(8, 6), dpi=100)
fig21.subplots_adjust(hspace=0.4)
ax211 = fig21.add_subplot(211)
ax212 = fig21.add_subplot(212)
fig22 = plt.figure(figsize=(8, 6), dpi=100)
fig22.subplots_adjust(hspace=0.4)
ax221 = fig22.add_subplot(211)
ax222 = fig22.add_subplot(212)
fig23 = plt.figure(figsize=(8, 6), dpi=100)
fig23.subplots_adjust(hspace=0.4)
ax231 = fig23.add_subplot(211)
ax232 = fig23.add_subplot(212)
fig24 = plt.figure(figsize=(8, 6), dpi=100)
fig24.subplots_adjust(hspace=0.4)
ax241 = fig24.add_subplot(211)
ax242 = fig24.add_subplot(212)

filepath = "./2/output_omega_q.csv"
df = pd.read_csv(filepath)

omega_x = df.iloc[:, 0]
omega_y = df.iloc[:, 1]
omega_z = df.iloc[:, 2]
q_0 = df.iloc[:, 3]
q_1 = df.iloc[:, 4]
q_2 = df.iloc[:, 5]
q_3 = df.iloc[:, 6]
omega_x_hat = df.iloc[:, 7]
omega_y_hat = df.iloc[:, 8]
omega_z_hat = df.iloc[:, 9]
q_0_hat = df.iloc[:, 10]
q_1_hat = df.iloc[:, 11]
q_2_hat = df.iloc[:, 12]
q_3_hat = df.iloc[:, 13]
P11 = df.iloc[:, 14]
P22 = df.iloc[:, 15]
P33 = df.iloc[:, 16]
P44 = df.iloc[:, 17]
P55 = df.iloc[:, 18]
P66 = df.iloc[:, 19]
P77 = df.iloc[:, 20]

dt = 0.01
dt_obs = 1
end_time = 50
t = np.arange(0, end_time + dt, dt)
t_len = len(t)

ax111.plot(t, omega_x, label="True", color="m", linestyle="-.", lw="1")
ax111.plot(t, omega_x_hat, label="Estimated", color="c")
ax112.plot(t, np.sqrt(P11), label="$\pm\sqrt{P_{11}}$", color="m", linestyle="-.", lw="1")
ax112.plot(t, -np.sqrt(P11), label="", color="m", linestyle="-.", lw="1")
ax112.plot(t, omega_x - omega_x_hat, label="$\Delta\omega_x$", color="c")
ax121.plot(t, omega_y, label="True", color="m", linestyle="-.", lw="1")
ax121.plot(t, omega_y_hat, label="Estimated", color="c")
ax122.plot(t, np.sqrt(P22), label="$\pm\sqrt{P_{22}}$", color="m", linestyle="-.", lw="1")
ax122.plot(t, -np.sqrt(P22), label="", color="m", linestyle="-.", lw="1")
ax122.plot(t, omega_y - omega_y_hat, label="$\Delta\omega_y$", color="c")
ax131.plot(t, omega_z, label="True", color="m", linestyle="-.", lw="1")
ax131.plot(t, omega_z_hat, label="Estimated", color="c")
ax132.plot(t, np.sqrt(P33), label="$\pm\sqrt{P_{33}}$", color="m", linestyle="-.", lw="1")
ax132.plot(t, -np.sqrt(P33), label="", color="m", linestyle="-.", lw="1")
ax132.plot(t, omega_z - omega_z_hat, label="$\Delta\omega_z$", color="c")

ax211.plot(t, q_0, label="True", color="m", linestyle="-.", lw="1")
ax211.plot(t, q_0_hat, label="Estimated", color="c")
ax212.plot(t, np.sqrt(P44), label="$\pm\sqrt{P_{44}}$", color="m", linestyle="-.", lw="1")
ax212.plot(t, -np.sqrt(P44), label="", color="m", linestyle="-.", lw="1")
ax212.plot(t, q_0 - q_0_hat, label="$\Delta q_0$", color="c")
ax221.plot(t, q_1, label="True", color="m", linestyle="-.", lw="1")
ax221.plot(t, q_1_hat, label="Estimated", color="c")
ax222.plot(t, np.sqrt(P55), label="$\pm\sqrt{P_{55}}$", color="m", linestyle="-.", lw="1")
ax222.plot(t, -np.sqrt(P55), label="", color="m", linestyle="-.", lw="1")
ax222.plot(t, q_1 - q_1_hat, label="$\Delta q_1$", color="c")
ax231.plot(t, q_2, label="True", color="m", linestyle="-.", lw="1")
ax231.plot(t, q_2_hat, label="Estimated", color="c")
ax232.plot(t, np.sqrt(P66), label="$\pm\sqrt{P_{66}}$", color="m", linestyle="-.", lw="1")
ax232.plot(t, -np.sqrt(P66), label="", color="m", linestyle="-.", lw="1")
ax232.plot(t, q_2 - q_2_hat, label="$\Delta q_2$", color="c")
ax241.plot(t, q_3, label="True", color="m", linestyle="-.", lw="1")
ax241.plot(t, q_3_hat, label="Estimated", color="c")
ax242.plot(t, np.sqrt(P77), label="$\pm\sqrt{P_{77}}$", color="m", linestyle="-.", lw="1")
ax242.plot(t, -np.sqrt(P77), label="", color="m", linestyle="-.", lw="1")
ax242.plot(t, q_3 - q_3_hat, label="$\Delta q_3$", color="c")


ax111.set_xlim(0, end_time)
ax111.set_ylim(-0.2, 0.2)
ax111.set_ylabel("$\omega_x\;\;$[rad/s]")
ax111.legend(loc="upper right").get_frame().set_alpha(1)
ax112.set_xlim(0, end_time)
ax112.set_ylim(-0.1, 0.1)
ax112.set_xlabel("$Time\;\;$[s]")
ax112.set_ylabel("$\Delta\omega_x\;\;$[rad/s]")
ax112.legend(loc="upper right").get_frame().set_alpha(1)
ax121.set_xlim(0, end_time)
ax121.set_ylim(1.875, 1.9)
ax121.set_ylabel("$\omega_y\;\;$[rad/s]")
ax121.legend(loc="upper right").get_frame().set_alpha(1)
ax122.set_xlim(0, end_time)
ax122.set_ylim(-0.1, 0.1)
ax122.set_xlabel("$Time\;\;$[s]")
ax122.set_ylabel("$\Delta\omega_y\;\;$[rad/s]")
ax122.legend(loc="upper right").get_frame().set_alpha(1)
ax131.set_xlim(0, end_time)
ax131.set_ylim(-0.2, 0.2)
ax131.set_ylabel("$\omega_z\;\;$[rad/s]")
ax131.legend(loc="upper right").get_frame().set_alpha(1)
ax132.set_xlim(0, end_time)
ax132.set_ylim(-0.1, 0.1)
ax132.set_xlabel("$Time\;\;$[s]")
ax132.set_ylabel("$\Delta\omega_z\;\;$[rad/s]")
ax132.legend(loc="upper right").get_frame().set_alpha(1)

ax211.set_xlim(0, end_time)
ax211.set_ylim(-1, 1)
ax211.set_ylabel("$q_0$")
ax211.legend(loc="upper right").get_frame().set_alpha(1)
ax212.set_xlim(0, end_time)
ax212.set_ylim(-0.1, 0.1)
ax212.set_xlabel("$Time\;\;$[s]")
ax212.set_ylabel("$\Delta q_0$")
ax212.legend(loc="upper right").get_frame().set_alpha(1)
ax221.set_xlim(0, end_time)
ax221.set_ylim(-0.2, 0.2)
ax221.set_ylabel("$q_1$")
ax221.legend(loc="upper right").get_frame().set_alpha(1)
ax222.set_xlim(0, end_time)
ax222.set_ylim(-0.1, 0.1)
ax222.set_xlabel("$Time\;\;$[s]")
ax222.set_ylabel("$\Delta q_1$")
ax222.legend(loc="upper right").get_frame().set_alpha(1)
ax231.set_xlim(0, end_time)
ax231.set_ylim(-1, 1)
ax231.set_ylabel("$q_2$")
ax231.legend(loc="upper right").get_frame().set_alpha(1)
ax232.set_xlim(0, end_time)
ax232.set_ylim(-0.1, 0.1)
ax232.set_xlabel("$Time\;\;$[s]")
ax232.set_ylabel("$\Delta q_2$")
ax232.legend(loc="upper right").get_frame().set_alpha(1)
ax241.set_xlim(0, end_time)
ax241.set_ylim(-0.2, 0.2)
ax241.set_ylabel("$q_3$")
ax241.legend(loc="upper right").get_frame().set_alpha(1)
ax242.set_xlim(0, end_time)
ax242.set_ylim(-0.1, 0.1)
ax242.set_xlabel("$Time\;\;$[s]")
ax242.set_ylabel("$\Delta q_3$")
ax242.legend(loc="upper right").get_frame().set_alpha(1)


plt.show()



# 推定誤差と√Pの比較（後者を１として規格化して図示した）
fig3 = plt.figure(figsize=(12, 8), dpi=100)
fig3.subplots_adjust(hspace=0.4)
ax31 = fig3.add_subplot(221)
ax32 = fig3.add_subplot(222)
ax33 = fig3.add_subplot(223)
fig4 = plt.figure(figsize=(12, 8), dpi=100)
fig4.subplots_adjust(hspace=0.4)
ax41 = fig4.add_subplot(221)
ax42 = fig4.add_subplot(222)
ax43 = fig4.add_subplot(223)
ax44 = fig4.add_subplot(224)

ax31.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{11}} / \sqrt{P_{11}}$", color="m", linestyle="-.", lw="1")
ax31.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax31.plot(t, (omega_x - omega_x_hat) / np.sqrt(P11), label="$\Delta\omega_x / \sqrt{P_{11}}$", color="c")
ax32.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{22}} / \sqrt{P_{22}}$", color="m", linestyle="-.", lw="1")
ax32.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax32.plot(t, (omega_y - omega_y_hat) / np.sqrt(P22), label="$\Delta\omega_y / \sqrt{P_{22}}$", color="c")
ax33.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{33}} / \sqrt{P_{33}}$", color="m", linestyle="-.", lw="1")
ax33.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax33.plot(t, (omega_z - omega_z_hat) / np.sqrt(P33), label="$\Delta\omega_z / \sqrt{P_{33}}$", color="c")
ax41.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{44}} / \sqrt{P_{44}}$", color="m", linestyle="-.", lw="1")
ax41.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax41.plot(t, (q_0 - q_0_hat) / np.sqrt(P11), label="$\Delta q_0 / \sqrt{P_{44}}$", color="c")
ax42.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{55}} / \sqrt{P_{55}}$", color="m", linestyle="-.", lw="1")
ax42.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax42.plot(t, (q_1 - q_1_hat) / np.sqrt(P22), label="$\Delta q_1 / \sqrt{P_{55}}$", color="c")
ax43.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{66}} / \sqrt{P_{66}}$", color="m", linestyle="-.", lw="1")
ax43.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax43.plot(t, (q_2 - q_2_hat) / np.sqrt(P33), label="$\Delta q_2 / \sqrt{P_{66}}$", color="c")
ax44.plot(t, [1 for i in range(t_len)], label="$\pm\sqrt{P_{77}} / \sqrt{P_{77}}$", color="m", linestyle="-.", lw="1")
ax44.plot(t, [-1 for i in range(t_len)], label="", color="m", linestyle="-.", lw="1")
ax44.plot(t, (q_3 - q_3_hat) / np.sqrt(P33), label="$\Delta q_3 / \sqrt{P_{77}}$", color="c")


ax31.set_xlim(0, end_time)
ax31.set_xlabel("$Time\;\;$[s]")
ax31.legend(loc="upper right").get_frame().set_alpha(1)
ax32.set_xlim(0, end_time)
ax32.set_xlabel("$Time\;\;$[s]")
ax32.legend(loc="lower right").get_frame().set_alpha(1)
ax33.set_xlim(0, end_time)
ax33.set_xlabel("$Time\;\;$[s]")
ax33.legend(loc="upper right").get_frame().set_alpha(1)
ax41.set_xlim(0, end_time)
ax41.set_xlabel("$Time\;\;$[s]")
ax41.legend(loc="upper right").get_frame().set_alpha(1)
ax42.set_xlim(0, end_time)
ax42.set_xlabel("$Time\;\;$[s]")
ax42.legend(loc="lower right").get_frame().set_alpha(1)
ax43.set_xlim(0, end_time)
ax43.set_xlabel("$Time\;\;$[s]")
ax43.legend(loc="upper right").get_frame().set_alpha(1)
ax44.set_xlim(0, end_time)
ax44.set_xlabel("$Time\;\;$[s]")
ax44.legend(loc="upper right").get_frame().set_alpha(1)

plt.show()


# 推定誤差が±√P内に収まるような状態変数の割合を時間軸でプロット
fig5 = plt.figure(figsize=(7, 4), dpi=100)
ax5 = fig5.add_subplot(111)

e_11 = [1 if abs(omega_x[i] - omega_x_hat[i]) / np.sqrt(P11[i]) < 1 else 0 for i in range(t_len)]
e_22 = [1 if abs(omega_y[i] - omega_y_hat[i]) / np.sqrt(P22[i]) < 1 else 0 for i in range(t_len)]
e_33 = [1 if abs(omega_z[i] - omega_z_hat[i]) / np.sqrt(P33[i]) < 1 else 0 for i in range(t_len)]
e_44 = [1 if abs(q_0[i] - q_0[i]) / np.sqrt(P44[i]) < 1 else 0 for i in range(t_len)]
e_55 = [1 if abs(q_1[i] - q_1[i]) / np.sqrt(P55[i]) < 1 else 0 for i in range(t_len)]
e_66 = [1 if abs(q_2[i] - q_2[i]) / np.sqrt(P66[i]) < 1 else 0 for i in range(t_len)]
e_77 = [1 if abs(q_3[i] - q_3[i]) / np.sqrt(P77[i]) < 1 else 0 for i in range(t_len)]

ax5.plot(t, [i.count(1) / 7 * 100 for i in zip(e_11, e_22, e_33, e_44, e_55, e_66, e_77)], color="c")
ax5.plot(t, [68.3 for i in range(t_len)], color="m", linestyle="-.", lw="1")
ax5.text(42, 69, "$1\sigma\simeq68.3\%$")
ax5.set_xlim(0, end_time)
ax5.set_xlabel("$Time\;\;$[s]")
ax5.set_ylabel("$\hat{x}\;\;in\;(-\sqrt{P_{ii}},\;\sqrt{P_{ii}})\;\;$[%]")

plt.show()

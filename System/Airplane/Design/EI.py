import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# y軸を指数表記にしたい
import matplotlib.ticker as ptick


filepath=r'C:\Users/tomot/Documents\航空宇宙\航空宇宙システム学製図\課題4\計算書.xlsx'
df=pd.read_excel(filepath, sheetname=4) #Pandas DataFrameに読込み

#データフレームからx軸とY軸データを読み込む
# #列で指定
# x=df[[0]]
# y=df[[1]]
# ラベルで指定する場合
x=df['y']
y=df['EI']

plt.rcParams["font.size"] = 14
fig = plt.figure()
ax = fig.add_subplot(111)


ax.plot(x,y,"-ok")
ax.set_xlim(0,5)
ax.set_ylim(0,2100000)
ax.set_xlabel("$y_w\;\;$[m]",fontsize=21)
ax.set_ylabel("$EI\;\;$[N m^2]",fontsize=21)

ax.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True)) 
ax.yaxis.offsetText.set_fontsize(14)
ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

plt.show()

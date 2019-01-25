import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


filepath=r'C:\Users/tomot/Documents\航空宇宙\航空宇宙システム学製図\課題4\計算書.xlsx'
df=pd.read_excel(filepath, sheetname=0) #Pandas DataFrameに読込み

#データフレームからx軸とY軸データを読み込む
# #列で指定
# x=df[[0]]
# y=df[[1]]
# ラベルで指定する場合
x=df['y']
y1=df['S']
y2=df['M']

plt.rcParams["font.size"] = 14

plt.plot(x,y1,"-ok")
plt.xlim(0,5)
plt.ylim(0,35000)
plt.xlabel("$y_w\;\;$[m]",fontsize=21)
plt.ylabel("$S\;\;$[N]",fontsize=21)
plt.show()

plt.plot(x,y2,"-ok")
plt.xlim(0,5)
plt.ylim(0,70000)
plt.xlabel("$y_w\;\;$[m]",fontsize=21)
plt.ylabel("$M\;\;$[N$\cdot$ m]",fontsize=21)
plt.show()

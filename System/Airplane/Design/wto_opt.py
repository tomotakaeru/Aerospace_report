"""
最大離陸重量WTOを，統計的関係式と，燃料重量・ペイロード重量で表される式
の二通りで導出し，相対誤差0.5%になるものを探す
"""

import numpy as np

wto_list = np.linspace(470000,480000,101)

for i in wto_list:
    if 0.9*pow(10, (np.log10(i)+0.163)/1.084 ) - (i*0.537-35614) <= 0.001*(i*0.537-35614):
        print(i)
print("end")

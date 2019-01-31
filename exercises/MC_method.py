import numpy as np

p_0 = np.array([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
T = np.array([
           [1/3,1/4,  0,1/4,  0,  0,  0,  0,  0],
           [1/3,1/4,1/3,  0,1/5,  0,  0,  0,  0],
           [  0,1/4,1/3,  0,  0,1/4,  0,  0,  0],
           [1/3,  0,  0,1/4,1/5,  0,1/3,  0,  0],
           [  0,1/4,  0,1/4,1/5,1/4,  0,1/4,  0],
           [  0,  0,1/3,  0,1/5,1/4,  0,  0,1/3],
           [  0,  0,  0,1/4,  0,  0,1/3,1/4,  0],
           [  0,  0,  0,  0,1/5,  0,1/3,1/4,1/3],
           [  0,  0,  0,  0,  0,1/4,  0,1/4,1/3]])


#1 繰り返し演算

def cycle_method(init,matrix,delta):
    p_n = init
    d_norm = np.linalg.norm(p_n) 
    while d_norm > delta:
        p_n0 = p_n
        p_n = np.dot(matrix,p_n0)
        d_norm = np.linalg.norm(p_n - p_n0)
    print(p_n)


#2 定常状態を求める
#行列Tの固有値1の固有ベクトルで各成分の和が1のものにpは収束する、これが定常値のpである

def steady_state_method(matrix):
    for index,la in enumerate(np.linalg.eig(matrix)[0]):
        if la==1:
            index_la1 = index
    V_la1 = np.linalg.eig(matrix)[1][:,index_la1]
    V_la1_nomalized = V_la1/sum(V_la1)
    print(V_la1_nomalized)


#3 モンテカルロ法

def mc_method(init,matrix,times):
    #初期位置
    for x,y in enumerate(init):
        if y==1:
            position_0 = x
    position = position_0

    #ある位置に何回いたかをカウントする配列
    position_cnt = np.zeros(len(init))
    position_cnt[position] += 1

    position_list = range(len(init))

    for i in range(times):
        weight = matrix[:,position]
        position = np.random.choice(position_list, p=weight)  
        if i>10000:
            position_cnt[position] += 1

    print(position_cnt/sum(position_cnt))



if __name__ == "__main__":
    cycle_method(p_0,T,0.000000001)
    steady_state_method(T)
    mc_method(p_0,T,1000000)    #pythonのforループのため20sほどかかる。小数点以下3桁まで他の2つの手法と一致。

"""
一般の連立一次方程式をGaussの消去法・Jacobiの反復法を用いて解く
入力：
    正則な係数行列 A
    定数ベクトル b
出力:
    解ベクトル x
"""
import numpy as np
import matplotlib.pyplot as plt

def direct_method(A,b): # Gaussの消去法
    A = A.copy()
    b = b.copy()
    N = A.shape[0]
    # 前進消去
    for i in range(N):
        for j in range(i+1,N):
            if A[j,i]==0:
                break
            else:
                w = A[j,i]/A[i,i]
                b[j] -= w*b[i]
                A[j,i] = 0
            for k in range(i+1,N):
                A[j,k] -= w*A[i,k]
    # 後退代入
    x = np.array([0. for i in range(N)])
    for i in reversed(range(N)):
        for j in range(i,N):
            b[i] -= A[i,j]*x[j]
            x[i] = b[i]/A[i,i]
    return x

def iteration_method(A,b,e=1e-3): # Jacobiの反復法(修正前)
    S = np.diag(A)
    A = A.copy()
    np.fill_diagonal(A,0)
    LU = A
    b = b.copy()
    N = A.shape[0]
    # x_next = np.array([b[i]/S[i] for i in range(N)])
    x_next = np.array([1.5,2.4,2.0,2.4,3.0])
    x = np.array([0. for i in range(N)])
    cnt=0
    S1=np.array([[2,0,0,0,0],[0,-8,0,0,0],[0,0,5,0,0],[0,0,0,3,0],[0,0,0,0,5]])
    eig = np.linalg.eig(np.dot(np.linalg.inv(S1),LU))[0]

    while np.linalg.norm(np.array(x_next)-np.array(x)) > e:
        x = x_next
        x_next = [(-np.dot(LU[i],x) + b[i])/S[i] for i in range(N)]
        cnt+=1
    return x_next,cnt,eig

def iteration_method_rev(A,b,e=1e-3): # Jacobiの反復法(修正後)
    S = np.diag(A)
    S = np.array(S)
    A = A.copy()
    np.fill_diagonal(A,0)
    LU = A
    b = b.copy()
    N = A.shape[0]
    x_next = np.array([b[i]/S[i] for i in range(N)])
    # x_next = np.array([1.5,2.4,2.0,2.4,3.0])
    x = np.array([0. for i in range(N)])
    cnt=0
    
    S[0]+=7
    S[1]+=0
    S[2]+=1
    S[3]+=3
    S[4]+=0
    LU[0,0]=-7
    LU[1,1]=0
    LU[2,2]=-1
    LU[3,3]=-3
    LU[4,4]=0

    S1=np.array([[16,0,0,0,0],[0,-8,0,0,0],[0,0,6,0,0],[0,0,0,6,0],[0,0,0,0,5]])
    eig = np.linalg.eig(np.dot(np.linalg.inv(S1),LU))[0]

    while np.linalg.norm(np.array(x_next)-np.array(x)) > e:
        x = x_next
        x_next = [(-np.dot(LU[i],x) + b[i])/S[i] for i in range(N)]
        cnt+=1
    return x_next,cnt,eig



if __name__ == "__main__":

    # 与えられた連立一次方程式
    A = np.array(
    [[2., 1,-1, 5,-1],
    [ 5,-8,-1,-1, 6],
    [-2, 1, 5, 4,-3],
    [ 1,-1,-4, 3, 2],
    [ 2, 1,-3,-1, 5]])
    b = np.array([12., 2,10, 4,12])

    x = direct_method(A,b)
    print("Gaussの消去法を用いた場合")
    print("解：",x)
    print("-------------------------------")

    print("Jacobiの反復法を用いた場合(修正前)")
    x,cnt,eig = iteration_method(A,b)
    print("解：",x)
    print("ステップ数：",cnt)
    print("固有値：",eig)
    print("固有値の絶対値の最大値：",max(abs(eig)))
    print("-------------------------------")

    print("Jacobiの反復法を用いた場合(修正後)")
    x,cnt,eig = iteration_method_rev(A,b)
    print("解：",x)
    print("ステップ数：",cnt)
    print("固有値：",eig)
    print("固有値の絶対値の最大値：",max(abs(eig)))

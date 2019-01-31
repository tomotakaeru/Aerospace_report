#剛体棒の二重振り子


from numpy import*
import matplotlib.pyplot as plt
import time


#系の状態
m=1.0
l=1.0
g=9.8
dt=0.1

#箱の用意
T = arange(0.0,20.0,dt)
N = T.shape[0]
X = zeros((N,4))

#初期条件[θ, θ', φ, φ']
X[0,:] =[pi/2, 0.0, 0.0, 0.0]


#加速度[θ'', φ'']を返す
def accel(list):
    x=list[0]
    x_dot=list[1]
    y=list[2]
    y_dot=list[3]

    right = -g/l * matrix([[3*sin(x),sin(y)]]).T +2*sin(x-y) * matrix([[-y_dot**2,x_dot**2]]).T
    A = array([[16/3, 2*cos(x-y)], [2*cos(x-y), 4/3]])
    if linalg.det(A) != 0:
        acc = dot(linalg.inv(A), right)
    else:
        print("逆行列が定義されません")
    return acc

#位置、速度を計算(4次のRunge-Kutta法)
def list_make(box):
    for n in range(N-1):
        a1 = accel(box[n,:])
        list1=[float(X[n,0] + X[n,1] * dt/2 + a1[0] * ((dt/2)**2)/2), float(X[n,1] + a1[0] * dt/2), 
        float(X[n,2] + X[n,3] * dt/2 + a1[1] * ((dt/2)**2)/2), float(X[n,3] + a1[1] * dt/2)]
        a2 = accel(list1)
        list2=[float(X[n,0] + X[n,1] * dt/2 + a2[0] * ((dt/2)**2)/2), float(X[n,1] + a2[0] * dt/2), 
        float(X[n,2] + X[n,3] * dt/2 + a2[1] * ((dt/2)**2)/2), float(X[n,3] + a2[1] * dt/2)]
        a3 = accel(list2)
        list3=[float(X[n,0] + X[n,1] * dt + a3[0] * (dt**2)/2), float(X[n,1] + a3[0] * dt), 
        float(X[n,2] + X[n,3] * dt + a3[1] * (dt**2)/2), float(X[n,3] + a3[1] * dt)]
        a4 = accel(list3)
        a=(a1+2*a2+2*a3+a4)/6
        
        X[n+1,1] = X[n,1] + a[0] * dt
        X[n+1,3] = X[n,3] + a[1] * dt
        X[n+1,0] = X[n,0] + X[n,1] * dt + a[0] * (dt**2)/2
        X[n+1,2] = X[n,2] + X[n,3] * dt + a[1] * (dt**2)/2

        if X[n+1,0] > pi:
            X[n+1,0]=-2*pi+X[n+1,0]
        if X[n+1,0] < -pi:
            X[n+1,0]=2*pi+X[n+1,0]
        if X[n+1,2] > pi:
            X[n+1,2]=-2*pi+X[n+1,2]
        if X[n+1,2] < -pi:
            X[n+1,2]=2*pi+X[n+1,2]

"""
#位置、速度を計算(Runge-Kutta法を用いない)
def list_make(box):
    for n in range(N-1):
        a=accel(box[n,:])
        X[n+1,1] = X[n,1] + a[0] * dt
        X[n+1,3] = X[n,3] + a[1] * dt
        X[n+1,0] = X[n,0] + X[n,1] * dt + a[0] * (dt**2)/2
        X[n+1,2] = X[n,2] + X[n,3] * dt + a[1] * (dt**2)/2

        if X[n+1,0] > pi:
            X[n+1,0]=-2*pi+X[n+1,0]
        if X[n+1,0] < -pi:
            X[n+1,0]=2*pi+X[n+1,0]
        if X[n+1,2] > pi:
            X[n+1,2]=-2*pi+X[n+1,2]
        if X[n+1,2] < -pi:
            X[n+1,2]=2*pi+X[n+1,2]
"""

#時系列に位置をプロット
def t_plot(box):
    list_make(box)
    plt.plot(T,X[:,0], label="θ")
    plt.plot(T,X[:,2], label="φ")
    plt.title("t - θ, φ plot")
    plt.ylim(-4,4)
    plt.legend()
    plt.draw()
    plt.pause(3)
    plt.close()

#位置座標をアニメーションプロット
def st_plot(box):
    list_make(box)
    x_s = [sin(X[n,0]) for n in range(N-1)]
    x_t = [-cos(X[n,0]) for n in range(N-1)]
    y_s = [sin(X[n,2]) for n in range(N-1)]
    y_t = [-cos(X[n,2]) for n in range(N-1)]
    y_s = [y_s[n] + x_s[n] for n in range(N-1)]
    y_t = [y_t[n] + x_t[n] for n in range(N-1)]

    rng_x=5 #θの表示幅
    rng_y=12 #φの表示幅
    mv=2 #刻み幅

    i=0
    theta=[0,0]
    fai=[0,0]
    while mv*i < rng_x:
        line, = plt.plot(x_s[:mv*i], x_t[:mv*i], linestyle="--", label="θ")
        line, = plt.plot(y_s[:mv*i], y_t[:mv*i], label="φ")
        theta[0] = x_s[mv*i]
        theta[1] = x_t[mv*i]
        fai[0] = y_s[mv*i]
        fai[1] = y_t[mv*i]
        line, = plt.plot([0,theta[0]], [0,theta[1]])
        line, = plt.plot([theta[0],fai[0]], [theta[1],fai[1]])
        plt.plot(x_s[mv*i], x_t[mv*i], marker="*")
        plt.plot(y_s[mv*i], y_t[mv*i], marker="*")

        plt.xlim(-2,2)
        plt.ylim(-2,2)
        plt.title("θ, φ animation plot")
        plt.legend()
        plt.draw()
        plt.pause(0.2)
        plt.clf()
        i+=1

    while mv*i < rng_y:
        line, = plt.plot(x_s[mv*i-rng_x:mv*i], x_t[mv*i-rng_x:mv*i], linestyle="--", label="θ")
        line, = plt.plot(y_s[:mv*i], y_t[:mv*i], label="φ")
        theta[0] = x_s[mv*i]
        theta[1] = x_t[mv*i]
        fai[0] = y_s[mv*i]
        fai[1] = y_t[mv*i]
        line, = plt.plot([0,theta[0]], [0,theta[1]])
        line, = plt.plot([theta[0],fai[0]], [theta[1],fai[1]])
        plt.plot(x_s[mv*i], x_t[mv*i], marker="*")
        plt.plot(y_s[mv*i], y_t[mv*i], marker="*")

        plt.xlim(-2,2)
        plt.ylim(-2,2)
        plt.title("θ, φ animation plot")
        plt.legend()
        plt.draw()
        plt.pause(0.2)
        plt.clf()
        i+=1

    for n in range(i,int(N/mv)-1):
        line, = plt.plot(x_s[mv*n-rng_x:mv*n], x_t[mv*n-rng_x:mv*n], linestyle="--", label="θ")
        line, = plt.plot(y_s[mv*n-rng_y:mv*n], y_t[mv*n-rng_y:mv*n], label="φ")
        theta[0] = x_s[mv*n]
        theta[1] = x_t[mv*n]
        fai[0] = y_s[mv*n]
        fai[1] = y_t[mv*n]
        line, = plt.plot([0,theta[0]], [0,theta[1]])
        line, = plt.plot([theta[0],fai[0]], [theta[1],fai[1]])
        plt.plot(x_s[mv*n], x_t[mv*n], marker="*")
        plt.plot(y_s[mv*n], y_t[mv*n], marker="*")

        plt.xlim(-2,2)
        plt.ylim(-2,2)
        plt.title("θ, φ animation plot")
        plt.legend()
        plt.draw()
        plt.pause(0.2)
        plt.clf()
    plt.close()



if __name__ == "__main__":
    t_plot(X)
    st_plot(X)
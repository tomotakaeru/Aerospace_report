"""
一般の二次元トラス構造について以下の入出力を返す，構造解析用モジュール
入力：
    node (全体座標系x座標[m], y座標[m], x変位境界条件:1, y変位境界条件:1)
    member (節点番号（数字小）, （数字大）, 断面積[cm2], ヤング率[GPa])
    load  (x荷重条件[kN], y荷重条件[kN])
出力:
    節点変位 d[cm]
    節点反力 f[kN]
    部材応力 sigma[MPa]
"""
import numpy as np
import matplotlib.pyplot as plt


class Node: # 節点の情報を持つ
    def __init__(self, node_data):
        self.x = node_data[0]
        self.y = node_data[1]
        self.x_condition = node_data[2]
        self.y_condition = node_data[3]

class Member: # 部材の情報を持つ
    def __init__(self, member_data):
        self.small_node = member_data[0]
        self.large_node = member_data[1]
        self.A = member_data[2]
        self.E = member_data[3]

class Load: # 節点荷重の情報を持つ
    def __init__(self, load_data):
        self.fx_condition = load_data[0]
        self.fy_condition = load_data[1]

class System: # 全体のシステム
    def __init__(self, node_init, member_init, load_init):
        self.node_num = len(node_init)
        self.member_num = len(member_init)
        self.load_num = len(load_init)
        self.nodes = [Node(node_init[i]) for i in range(self.node_num)]
        self.members = [Member(member_init[i]) for i in range(self.member_num)]
        self.loads = [Load(load_init[i]) for i in range(self.load_num)]

        self.small_node_index = [int(self.members[i].small_node)-1 for i in range(self.member_num)]
        self.large_node_index = [int(self.members[i].large_node)-1 for i in range(self.member_num)]
        self.L = [ np.sqrt( (self.nodes[self.large_node_index[i]].y - self.nodes[self.small_node_index[i]].y)**2
         + (self.nodes[self.large_node_index[i]].x - self.nodes[self.small_node_index[i]].x)**2 ) for i in range(self.member_num) ]
        self.cos = [ (self.nodes[self.large_node_index[i]].x - self.nodes[self.small_node_index[i]].x) / self.L[i] for i in range(self.member_num) ]
        self.sin = [ (self.nodes[self.large_node_index[i]].y - self.nodes[self.small_node_index[i]].y) / self.L[i] for i in range(self.member_num) ]
        self.stage = 0

    def calc_part_stiffness_matrix(self): # 各部材の剛性マトリックスを求める
        self.stage+=1
        
        self.K_p = [ self.members[i].E * self.members[i].A / self.L[i]
        * np.array([[self.cos[i]**2, self.cos[i]*self.sin[i], -self.cos[i]**2, -self.cos[i]*self.sin[i]],
        [self.cos[i]*self.sin[i], self.sin[i]**2, -self.cos[i]*self.sin[i], -self.sin[i]**2],
        [-self.cos[i]**2, -self.cos[i]*self.sin[i], self.cos[i]**2, self.cos[i]*self.sin[i]],
        [-self.cos[i]*self.sin[i], -self.sin[i]**2, self.cos[i]*self.sin[i], self.sin[i]**2]
        ]) for i in range(self.member_num) ]

    def calc_overall_stiffness_matrix(self): # 全体剛性マトリックスを求める
        if self.stage<1:
            print("WARNING: skip some necessary process")
        else:
            self.stage+=1
            self.K_a = np.zeros((self.node_num*2,self.node_num*2))
            for i in range(self.member_num):
                self.K_a[self.small_node_index[i]*2:self.small_node_index[i]*2+2, self.small_node_index[i]*2:self.small_node_index[i]*2+2] += self.K_p[i][0:2,0:2]
                self.K_a[self.small_node_index[i]*2:self.small_node_index[i]*2+2, self.large_node_index[i]*2:self.large_node_index[i]*2+2] += self.K_p[i][0:2,2:4]
                self.K_a[self.large_node_index[i]*2:self.large_node_index[i]*2+2, self.small_node_index[i]*2:self.small_node_index[i]*2+2] += self.K_p[i][2:4,0:2]
                self.K_a[self.large_node_index[i]*2:self.large_node_index[i]*2+2, self.large_node_index[i]*2:self.large_node_index[i]*2+2] += self.K_p[i][2:4,2:4]

    def solve_equation(self):
        if self.stage<2:
            print("WARNING: skip some necessary process")
        else:
            self.stage+=1
            cnt=0

            # 変位境界条件があるnodeのindexを作る
            condition_index=[]
            for i in range(self.node_num):
                if self.nodes[i].x_condition != 0:
                    condition_index.append(cnt)
                cnt+=1
                if self.nodes[i].y_condition != 0:
                    condition_index.append(cnt)
                cnt+=1
            
            # 分割された行列を得る
            self.K_a_origin = self.K_a.copy()
            self.K_a = np.vstack((self.K_a, self.K_a[condition_index,:]))
            self.K_a = np.delete(self.K_a, condition_index, 0)
            self.K_a = np.hstack((self.K_a, self.K_a[:,condition_index]))
            self.K_a = np.delete(self.K_a, condition_index, 1)
            Kff = self.K_a[0:self.node_num*2-len(condition_index), 0:self.node_num*2-len(condition_index)]
            Ksf = self.K_a[0:self.node_num*2-len(condition_index), self.node_num*2-len(condition_index):self.node_num*2]
            Kfs = self.K_a[self.node_num*2-len(condition_index):self.node_num*2, 0:self.node_num*2-len(condition_index)]
            Kss = self.K_a[self.node_num*2-len(condition_index):self.node_num*2, self.node_num*2-len(condition_index):self.node_num*2]
            ds = np.c_[np.zeros(len(condition_index))]
            self.f = np.c_[np.hstack(tuple([np.array([self.loads[i].fx_condition, self.loads[i].fy_condition]) for i in range(self.load_num)]))]
            ff = np.delete(self.f, condition_index, 0)

            # 剛性方程式を解く
            self.df = np.linalg.inv(Kff).dot(ff-Ksf.dot(ds))
            self.fs = Kfs.dot(self.df) + Kss.dot(ds)

            # 境界条件の値も加えて，d,fベクトルとする
            self.d = self.df.copy()
            for i in range(len(condition_index)):
                self.d = np.insert(self.d, condition_index[i], 0, 0)
                self.f[condition_index[i]] = self.fs[i]

            # 計算誤差
            self.error = self.K_a_origin.dot(self.d) - self.f

    def calc_stress(self): # 各部材の応力を求める
        if self.stage<3:
            print("WARNING: skip some necessary process")
        else:
            self.stage+=1
            self.sigma = [ 10* ( self.members[i].E / self.L[i]
             * np.array([-self.cos[i], -self.sin[i], self.cos[i], self.sin[i]]).dot(np.array([ self.d[self.small_node_index[i]*2], self.d[self.small_node_index[i]*2+1], self.d[self.large_node_index[i]*2], self.d[self.large_node_index[i]*2+1] ])) )[0]
             for i in range(self.member_num) ]

    def print_result(self):
        print("全体剛性マトリックス：")
        print(self.K_a_origin)
        print("-----------------------------------------------------------------------------------")
        print("各節点変位[cm] (x方向, y方向)：")
        for i in range(self.node_num):
            print("節点", i+1, "\t", (self.d[i*2][0],self.d[i*2+1][0]))
        print("-----------------------------------------------------------------------------------")
        print("各節点荷重[kN] (x方向, y方向)：")
        for i in range(self.node_num):
            print("節点", i+1, "\t", (self.f[i*2][0],self.f[i*2+1][0]))
        print("-----------------------------------------------------------------------------------")
        print("計算誤差：")
        print(self.error.T)
        print("-----------------------------------------------------------------------------------")
        print("各部材応力[MPa]：")
        for i in range(self.member_num):
            print("部材", i+1, "\t", self.sigma[i])

    def plot_result(self, displace_scale=1.0, load_scale=1.0): # displace_scaleで変位，load_scaleで荷重矢印を拡大する
        # 変形前の形状
        for i in range(self.node_num):
            if self.nodes[i].x_condition != 0:
                plt.plot(self.nodes[i].x, self.nodes[i].y, "sk", markersize=9)
            else:
                if self.nodes[i].y_condition != 0:
                    plt.plot(self.nodes[i].x, self.nodes[i].y, "^k", markersize=9)
                else:
                    plt.plot(self.nodes[i].x, self.nodes[i].y, "ok")
        for i in range(self.member_num):
            p1, = plt.plot([self.nodes[self.small_node_index[i]].x, self.nodes[self.large_node_index[i]].x] , [self.nodes[self.small_node_index[i]].y, self.nodes[self.large_node_index[i]].y], "k")
        
        # 変形後の形状
        for i in range(self.node_num):
            if self.nodes[i].x_condition != 0:
                plt.plot(self.nodes[i].x+self.d[i*2]/100*displace_scale, self.nodes[i].y+self.d[i*2+1][0]/100*displace_scale, "sr", markersize=9)
            else:
                if self.nodes[i].y_condition != 0:
                    plt.plot(self.nodes[i].x+self.d[i*2]/100*displace_scale, self.nodes[i].y+self.d[i*2+1][0]/100*displace_scale, "^r", markersize=9)
                else:
                    plt.plot(self.nodes[i].x+self.d[i*2]/100*displace_scale, self.nodes[i].y+self.d[i*2+1][0]/100*displace_scale, "or")
        for i in range(self.member_num):
            p2, = plt.plot( [self.nodes[self.small_node_index[i]].x+self.d[self.small_node_index[i]*2]/100*displace_scale , self.nodes[self.large_node_index[i]].x+self.d[self.large_node_index[i]*2]/100*displace_scale] ,
             [self.nodes[self.small_node_index[i]].y+self.d[self.small_node_index[i]*2+1]/100*displace_scale , self.nodes[self.large_node_index[i]].y+self.d[self.large_node_index[i]*2+1]/100*displace_scale], "r")

        # 荷重を矢印で表現
        for i in range(self.node_num):
            plt.quiver(self.nodes[i].x+self.d[i*2]/100*displace_scale, self.nodes[i].y+self.d[i*2+1][0]/100*displace_scale,
             self.f[i*2]/200*load_scale, self.f[i*2+1]/200*load_scale, angles='xy',scale_units='xy',scale=1)

        # プロットの設定
        plt.axes().set_aspect('equal', 'datalim')
        plt.legend([p1, p2], ["original", "result"], prop={'size':15,})
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.show()


def main(node_data, member_data, load_data):
    system = System(node_data, member_data, load_data)
    system.calc_part_stiffness_matrix()
    system.calc_overall_stiffness_matrix()
    system.solve_equation()
    system.calc_stress()
    system.print_result()
    system.plot_result(3)


if __name__ == "__main__":

    # 与えられた条件
    node_init = np.array( # (全体座標系x座標[m], y座標[m], x変位境界条件:1, y変位境界条件:1)
    [[0., 0, 1, 1],
    [10, 10, 0, 0],
    [10,  0, 0, 0],
    [20, 10, 0, 0],
    [20,  0, 0, 0],
    [30,  0, 0, 1]])

    member_init = np.array( # (節点番号（数字小）, （数字大）, 断面積[cm2], ヤング率[GPa])
    [[1., 2, 14.142, 70],
    [ 1,  3,     10, 70],
    [ 2,  3,     10, 70],
    [ 2,  4,     10, 70],
    [ 3,  5,     10, 70],
    [ 3,  4, 14.142, 70],
    [ 2,  5, 14.142, 70],
    [ 4,  5,     10, 70],
    [ 5,  6,     10, 70],
    [ 4,  6, 14.142, 70]])

    load_init = np.array( # (x荷重条件[kN], y荷重条件[kN])
    [[0.,   0],
    [ 0,    0],
    [ 0, -900],
    [ 0,    0],
    [ 0, -500],
    [ 0,    0]])


    main(node_init, member_init, load_init)

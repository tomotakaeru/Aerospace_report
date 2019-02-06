#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>


/* 計算条件 */
const double Re = 70.0;  // Reynolds Number
const double cfl = 0.2;  // CFL Number
// SOR Pamameters
const double omega = 1.00;
const int max_itr = 100;
const double error = 0.0001;
// No. of Time Steps
const int nlast = 10000;  // 5min/1000

/* 計算格子の設定 */
// x-grid parameters
const int mx = 401;   // x軸格子点数(1~401),x=30.0がmx=401に対応
const int i_1 = 96;   // x軸格子点数基準での,角柱の左端
const int i_2 = 106;  // x軸格子点数基準での,角柱の右端
// y-grid parameters
const int my = 201;   // y軸格子点数(1~201)
const int j_1 = 96;   // y軸格子点数基準での,角柱の下端
const int j_2 = 106;  // y軸格子点数基準での,角柱の上端
// set delta x,y,t
const double dx = 1.0 / (i_2 - i_1);
const double dy = 1.0 / (j_2 - j_1);
const double dt = cfl * fmin(dx, dy);
// xy-grid system
double x[mx + 1];	// はじめに一つ余分に
double y[my + 1];	// はじめに一つ余分に
double u[mx + 2][my + 2], v[mx + 2][my + 2], p[mx + 2][my + 2];	// 外側に一つずつ余分に
double cp[mx + 2][my + 2];	// 外側に一つずつ余分に

/* 初期条件として一様流を設定 */
void set_init_condition(double u[mx + 2][my + 2], double v[mx + 2][my + 2], double p[mx + 2][my + 2]) {
    for (int i = 1; i <= mx; i++) {
        for (int j = 1; j <= my; j++) {
            u[i][j] = 1.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }
}

/* 圧力場の境界条件 */
void set_boundary_condition_p(double p[mx + 2][my + 2]) {
	// 流入出条件
    for (int j = 1; j <= my; j++) {
        p[1][j] = 0.0;
        p[mx][j] = 0.0;
    }
    for (int i = 1; i <= mx; i++) {
        p[i][1] = 0.0;
        p[i][my] = 0.0;
    }
	// 角柱周りの条件（四隅＆四辺）    
    p[i_1][j_1] = p[i_1 - 1][j_1 - 1];
    p[i_1][j_2] = p[i_1 - 1][j_2 + 1];
    p[i_2][j_1] = p[i_2 + 1][j_1 - 1];
    p[i_2][j_2] = p[i_2 + 1][j_2 + 1];
    for (int j = j_1 + 1; j <= j_2 - 1; j++) {
        p[i_1][j] = p[i_1 - 1][j];
        p[i_2][j] = p[i_2 + 1][j];
    }
    for (int i = i_1 + 1; i <= i_2 - 1; i++) {
        p[i][j_1] = p[i][j_1 - 1];
        p[i][j_2] = p[i][j_2 + 1];
    }
}

/* 速度場の境界条件 */
void set_boundary_condition_V(double u[mx + 2][my + 2], double v[mx + 2][my + 2]) {
	// 流入出条件
    for (int j = 1; j <= my; j++) {
        u[1][j] = 1.0;
        v[1][j] = 0.0;
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        u[mx][j] = 2.0 * u[mx - 1][j] - u[mx - 2][j];
        v[mx][j] = 2.0 * v[mx - 1][j] - v[mx - 2][j];
        u[mx + 1][j] = 2.0 * u[mx][j] - u[mx - 1][j];
        v[mx + 1][j] = 2.0 * v[mx][j] - v[mx - 1][j];
    }
    for (int i = 1; i <= mx; i++) {
        u[i][1] = 2.0 * u[i][2] - u[i][3];
        v[i][1] = 2.0 * v[i][2] - v[i][3];
        u[i][0] = 2.0 * u[i][1] - u[i][2];
        v[i][0] = 2.0 * v[i][1] - v[i][2];
        u[i][my] = 2.0 * u[i][my - 1] - u[i][my - 2];
        v[i][my] = 2.0 * v[i][my - 1] - v[i][my - 2];
        u[i][my + 1] = 2.0 * u[i][my] - u[i][my - 1];
        v[i][my + 1] = 2.0 * v[i][my] - v[i][my - 1];
    }
	// 角柱周りの条件（四隅＆四辺）
    for (int i = i_1; i <= i_2; i++) {
        for (int j = j_1; j <= j_2; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

/* 圧力場を求める（ポワソン方程式をSORで解く） */
void solve_poisson_eq(double u[mx + 2][my + 2], double v[mx + 2][my + 2], double p[mx + 2][my + 2], double dx, double dy, double dt) {
    double rhs[mx + 2][my + 2];	// 方程式の右辺
    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {	// 角柱内は計算を飛ばす
                continue;
            }
            double ux = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
            double uy = (u[i][j + 1] - u[i][j - 1]) / (2.0 * dy);
            double vx = (v[i + 1][j] - v[i - 1][j]) / (2.0 * dx);
            double vy = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy);
            rhs[i][j] = (ux + vy) / dt - (ux * ux + 2.0 * uy * vx + vy * vy);
        }
    }

    for (int itr = 1; itr <= max_itr; itr++) {
        double res = 0.0;
        for (int i = 2; i <= mx - 1; i++) {
            for (int j = 2; j <= my - 1; j++) {
                if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                    continue;
                }
                double dp = (p[i + 1][j] + p[i - 1][j]) / (dx * dx) + (p[i][j + 1] + p[i][j - 1]) / (dy * dy) - rhs[i][j];
                dp = dp / (2.0 / (dx * dx) + 2.0 / (dy * dy)) - p[i][j];
                res = res + dp * dp;
                p[i][j] = p[i][j] + omega * dp;
            }
        }
        set_boundary_condition_p(p);
        res = sqrt(res / double(mx * my));
        if (res < error) break;	//残差が十分小さくなる(収束する)までループ
    }
}

/* 速度場を求める（運動量式をKawamura-Kuwaharaスキームなどで解く） */
void solve_velocity_eq(double u[mx + 2][my + 2], double v[mx + 2][my + 2], double p[mx + 2][my + 2], double dx, double dy, double dt) {
    double urhs[mx + 2][my + 2], vrhs[mx + 2][my + 2];
	// 圧力勾配項
    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                continue;
            }
            urhs[i][j] = -(p[i + 1][j] - p[i - 1][j]) / (2.0 * dx);
            vrhs[i][j] = -(p[i][j + 1] - p[i][j - 1]) / (2.0 * dy);
        }
    }
	// 粘性項
    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                continue;
            }
            urhs[i][j] += (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / (Re * dx * dx) + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * dy * dy);
            vrhs[i][j] += (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / (Re * dx * dx) + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / (Re * dy * dy);
        }
    }

	// 移流項
    for (int j = j_1 + 1; j <= j_2 - 1; j++) {	// 角柱内部まで外挿
        u[i_1 + 1][j] = 2.0 * u[i_1][j] - u[i_1 - 1][j];
        u[i_2 - 1][j] = 2.0 * u[i_2][j] - u[i_2 + 1][j];
        v[i_1 + 1][j] = 2.0 * v[i_1][j] - v[i_1 - 1][j];
        v[i_2 - 1][j] = 2.0 * v[i_2][j] - v[i_2 + 1][j];
    }

    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                continue;
            }
            urhs[i][j] += - u[i][j] * (-u[i + 2][j] + 8.0 * (u[i + 1][j] - u[i - 1][j]) + u[i - 2][j]) / (12.0 * dx)
                - abs(u[i][j]) * (u[i + 2][j] - 4.0 * u[i + 1][j] + 6.0 * u[i][j] - 4.0 * u[i - 1][j] + u[i - 2][j]) / (4.0 * dx);
            vrhs[i][j] += - u[i][j] * (-v[i + 2][j] + 8.0 * (v[i + 1][j] - v[i - 1][j]) + v[i - 2][j]) / (12.0 * dx)
                - abs(u[i][j]) * (v[i + 2][j] - 4.0 * v[i + 1][j] + 6.0 * v[i][j] - 4.0 * v[i - 1][j] + v[i - 2][j]) / (4.0 * dx);
        }
    }

    for (int i = i_1 + 1; i <= i_2 - 1; i++) {	// 角柱内部まで外挿
        u[i][j_1 + 1] = 2.0 * u[i][j_1] - u[i][j_1 - 1];
        u[i][j_2 - 1] = 2.0 * u[i][j_2] - u[i][j_2 + 1];
        v[i][j_1 + 1] = 2.0 * v[i][j_1] - v[i][j_1 - 1];
        v[i][j_2 - 1] = 2.0 * v[i][j_2] - v[i][j_2 + 1];
    }
    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                continue;
            }
            urhs[i][j] += - v[i][j] * (-u[i][j + 2] + 8.0 * (u[i][j + 1] - u[i][j - 1]) + u[i][j - 2]) / (12.0 * dy)
                - abs(v[i][j]) * (u[i][j + 2] - 4.0 * u[i][j + 1] + 6.0 * u[i][j] - 4.0 * u[i][j - 1] + u[i][j - 2]) / (4.0 * dy);
            vrhs[i][j] += - v[i][j] * (-v[i][j + 2] + 8.0 * (v[i][j + 1] - v[i][j - 1]) + v[i][j - 2]) / (12.0 * dy)
                - abs(v[i][j]) * (v[i][j + 2] - 4.0 * v[i][j + 1] + 6.0 * v[i][j] - 4.0 * v[i][j - 1] + v[i][j - 2]) / (4.0 * dy);
        }
    }

	// 更新
    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if (i_1 < i && i < i_2 && j_1 < j && j < j_2) {
                continue;
            }
            u[i][j] += dt * urhs[i][j];
            v[i][j] += dt * vrhs[i][j];
        }
    }
}

/* csv出力 */
int file_write_p(double p[mx + 2][my + 2], const char *file_name) {
    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = my; i >= 1; i--) {
	for (int j = 1; j <= mx; j++) {
	    fprintf(fp, "%f", p[j][i]);
	    if (j == mx) {
		fprintf(fp, "\n");
	    } else {
		fprintf(fp, ",");
	    }
	}
    }
    fclose(fp);
    return 0;
}
int file_write_cdclcp(double cd, double cl, double cp1, double cp2) {
    FILE *fp;
    fp = fopen("./rectangular/data/cdclcp.csv", "a");
    fprintf(fp, "%f,%f,%f,%f\n", cd, cl, cp1, cp2);
    fclose(fp);
    return 0;
}



int main() {
    set_init_condition(u, v, p);
    set_boundary_condition_p(p);
    set_boundary_condition_V(u, v);

    std::chrono::system_clock::time_point start_time, end_time;  // 時間計測用変数を用意
    start_time = std::chrono::system_clock::now();  // 計測始まり

    for (int n = 1; n <= nlast; n++) {
		solve_poisson_eq(u, v, p, dx, dy, dt);
		set_boundary_condition_p(p);
		solve_velocity_eq(u, v, p, dx, dy, dt);
		set_boundary_condition_V(u, v);

		double cd = 0.0;
		double cl = 0.0;
		for (int j = j_1; j <= j_2 - 1; j++) {
			double cpfore = (2 * p[i_1][j] + 2 * p[i_1][j + 1]) / 2;
			double cpback = (2 * p[i_2][j] + 2 * p[i_2][j + 1]) / 2;
			cd += (cpfore - cpback) * dy;
		}
		for (int i = i_1; i <= i_2 - 1; i++) {
			double cpbtm = (2 * p[i][j_1] + 2 * p[i + 1][j_1]) / 2;
			double cptop = (2 * p[i][j_2] + 2 * p[i + 1][j_2]) / 2;
			cl += (cpbtm - cptop) * dx;
		}
		for (int i = 1; i <= mx; i++) {
			for (int j = 1; j <= my; j++) {
				cp[i][j] = 2 * p[i][j];
			}
		}
		double cp1 = 2 * p[i_2 + i_2 - i_1][j_1];
		double cp2 = 2 * p[i_2 + i_2 - i_1][j_2];

		end_time = std::chrono::system_clock::now();  // 計測終わり
		double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time -start_time).count() / 1000.0;	 // かかった時間[s]
		std::cout << "step:" << n << "/" << nlast << " elapsed_time:" << elapsed_time << "[s]" << std::endl;

		std::string s = "./rectangular/data/p" + std::to_string(n) + ".csv";
		const char *cs = s.data();
		file_write_p(cp, cs);
		file_write_cdclcp(cd, cl, cp1, cp2);
    }

    return 0;
}

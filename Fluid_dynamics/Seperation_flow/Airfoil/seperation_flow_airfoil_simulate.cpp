#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


/* 計算条件 */
const double Re = 70.0; // Reynolds Number
const double cfl = 0.2; // CFL Number
// SOR Pamameters
const double omega = 1.00;
const int max_itr = 100;
const double error = 0.001;
// No. of Time Steps
const int nlast = 10000;    // 1min / 100step

/* 計算格子の設定 */
// x-grid parameters
const int mx = 801; // x軸格子点数(1~801)
const int i_1 = 171;    // x軸格子点数基準での，翼の大体の左端
const int i_2 = 231;    // x軸格子点数基準での，翼の大体の右端
// y-grid parameters
const int my = 401; // y軸格子点数(1~401)
const int j_1 = 198;    // y軸格子点数基準での，翼の大体の下端
const int j_2 = 216;    // y軸格子点数基準での，翼の大体の上端
// set delta x,y,t
const double dx = 1.0 / (j_2 - j_1);
const double dy = 1.0 / (j_2 - j_1);
const double dt = cfl * fmin(dx, dy);
// xy-grid system
double u[mx + 2][my + 2], v[mx + 2][my + 2], p[mx + 2][my + 2];	// 外側に一つずつ余分に
double cp[mx + 2][my + 2];	// 外側に一つずつ余分に

/* vectorからvalueを探し，indexを返す．無かったら-1． */
int find_index(int value, std::vector<int> vec) {
    auto itr = std::find(vec.begin(), vec.end(), value);
    size_t index = std::distance(vec.begin(), itr);
    if(index == vec.size()) {
        return -1;
    }
    return index;
}

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
void set_boundary_condition_p(double p[mx + 2][my + 2], std::vector<int> wing_wall_x, std::vector<int> wing_wall_y, std::vector<int> wing_x, std::vector<int> wing_y) {
	// 流入出条件
    for (int j = 1; j <= my; j++) {
        p[1][j] = 0.0;
        p[mx][j] = 0.0;
    }
    for (int i = 1; i <= mx; i++) {
        p[i][1] = 0.0;
        p[i][my] = 0.0;
    }
    // 翼内部の条件
    for( const auto& el : wing_x ){
        size_t index = &el - &wing_x[0];
        p[wing_x[index]][wing_y[index]] = 0.0;
    }
	// 翼周りの条件
    for( const auto& el : wing_wall_x ){
        size_t index = &el - &wing_wall_x[0];
        if (index<72) {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]][wing_wall_y[index]+1];
        } else if (index==72) {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]-1][wing_wall_y[index]];
        } else {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]][wing_wall_y[index]-1];
        }
    }
    for( const auto& el : wing_wall_x ){    // 上下方向に連続する二点の場合に内側の点が外の値を参照するべく二回回す
        size_t index = &el - &wing_wall_x[0];
        if (index<72) {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]][wing_wall_y[index]+1];
        } else if (index==72) {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]-1][wing_wall_y[index]];
        } else {
            p[wing_wall_x[index]][wing_wall_y[index]] = p[wing_wall_x[index]][wing_wall_y[index]-1];
        }
    }
}

/* 速度場の境界条件 */
void set_boundary_condition_V(double u[mx + 2][my + 2], double v[mx + 2][my + 2], std::vector<int> wing_x, std::vector<int> wing_y) {
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
    // 翼内部の条件
    for( const auto& el : wing_x ){
        size_t index = &el - &wing_x[0];
        u[wing_x[index]][wing_y[index]] = 0.0;
        v[wing_x[index]][wing_y[index]] = 0.0;
    }
}

/* 圧力場を求める（ポワソン方程式を緩和法で解く） */
void solve_poisson_eq(double u[mx + 2][my + 2], double v[mx + 2][my + 2], double p[mx + 2][my + 2], double dx, double dy, double dt, std::vector<int> wing_wall_x, std::vector<int> wing_wall_y, std::vector<int> wing_x, std::vector<int> wing_y) {
    double rhs[mx + 2][my + 2];	// 方程式の右辺をまず計算
    for (int i = 2; i <= mx - 1; i++) {
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
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
            int index = find_index(i, wing_x);
            for (int j = 2; j <= my - 1; j++) {
                if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
                    continue;
                }
                double dp = (p[i + 1][j] + p[i - 1][j]) / (dx * dx) + (p[i][j + 1] + p[i][j - 1]) / (dy * dy) - rhs[i][j];
                dp = dp / (2.0 / (dx * dx) + 2.0 / (dy * dy)) - p[i][j];
                res += dp * dp;
                p[i][j] += omega * dp;
            }
        }
        set_boundary_condition_p(p, wing_wall_x, wing_wall_y, wing_x, wing_y);  // 境界条件を適用
        res = sqrt(res / double(mx * my));
        if (res < error) break;	// 残差が十分小さくなる(収束する)までループ
    }
}

/* 速度場を求める（運動量式をKawamura-Kuwaharaスキームなどで解く） */
void solve_velocity_eq(double u[mx + 2][my + 2], double v[mx + 2][my + 2], double p[mx + 2][my + 2], double dx, double dy, double dt, std::vector<int> wing_wall_x, std::vector<int> wing_wall_y, std::vector<int> wing_x, std::vector<int> wing_y) {
    double urhs[mx + 2][my + 2], vrhs[mx + 2][my + 2];
	// 圧力勾配項
    for (int i = 2; i <= mx - 1; i++) {
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
                continue;
            }
            urhs[i][j] = -(p[i + 1][j] - p[i - 1][j]) / (2.0 * dx);
            vrhs[i][j] = -(p[i][j + 1] - p[i][j - 1]) / (2.0 * dy);
        }
    }
	// 粘性項
    for (int i = 2; i <= mx - 1; i++) {
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
                continue;
            }
            urhs[i][j] += (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / (Re * dx * dx) + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * dy * dy);
            vrhs[i][j] += (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / (Re * dx * dx) + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / (Re * dy * dy);
        }
    }

	// 移流項
    // x方向の差分
    for (const auto& i: wing_wall_x) {	// 翼内部まで外挿
        int index = find_index(i, wing_wall_x);

        if (find_index(i+1, wing_x) != -1 && find_index(i+1, wing_wall_x) == -1) {  // 一つ右の格子点が翼内部かつ翼壁でないとき
            int j = wing_wall_y[index];
            u[i+1][j] = 2.0 * u[i][j] - u[i-1][j];
            v[i+1][j] = 2.0 * v[i][j] - v[i-1][j];
        }
        if (find_index(i-1, wing_x) != -1 && find_index(i-1, wing_wall_x) == -1) {  // 一つ左の格子点が翼内部かつ翼壁でないとき
            int j = wing_wall_y[index];
            u[i-1][j] = 2.0 * u[i][j] - u[i+1][j];
            v[i-1][j] = 2.0 * v[i][j] - v[i+1][j];
        }
    }
    for (int i = 2; i <= mx - 1; i++) {
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
                continue;
            }
            urhs[i][j] += - u[i][j] * (-u[i + 2][j] + 8.0 * (u[i + 1][j] - u[i - 1][j]) + u[i - 2][j]) / (12.0 * dx)
                - abs(u[i][j]) * (u[i + 2][j] - 4.0 * u[i + 1][j] + 6.0 * u[i][j] - 4.0 * u[i - 1][j] + u[i - 2][j]) / (4.0 * dx);
            vrhs[i][j] += - u[i][j] * (-v[i + 2][j] + 8.0 * (v[i + 1][j] - v[i - 1][j]) + v[i - 2][j]) / (12.0 * dx)
                - abs(u[i][j]) * (v[i + 2][j] - 4.0 * v[i + 1][j] + 6.0 * v[i][j] - 4.0 * v[i - 1][j] + v[i - 2][j]) / (4.0 * dx);
        }
    }

    // y方向の差分
    for (const auto& j: wing_wall_y) {	// 翼内部まで外挿
        int index = find_index(j, wing_wall_y);

        if (find_index(j+1, wing_y) != -1 && find_index(j+1, wing_wall_y) == -1) {  // 一つ上の格子点が翼内部かつ翼壁でないとき
            int i = wing_wall_x[index];
            u[i][j+1] = 2.0 * u[i][j] - u[i][j-1];
            v[i][j+1] = 2.0 * v[i][j] - v[i][j-1];
        }
        if (find_index(j-1, wing_y) != -1 && find_index(j-1, wing_wall_y) == -1) {  // 一つ下の格子点が翼内部かつ翼壁でないとき
            int i = wing_wall_x[index];
            u[i][j-1] = 2.0 * u[i][j] - u[i][j+1];
            v[i][j-1] = 2.0 * v[i][j] - v[i][j+1];
        }
    }
    for (int i = 2; i <= mx - 1; i++) {
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
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
        int index = find_index(i, wing_x);
        for (int j = 2; j <= my - 1; j++) {
            if (index != -1 && j == wing_y[index]) {	// 翼内部は計算を飛ばす
                continue;
            }
            u[i][j] += dt * urhs[i][j];
            v[i][j] += dt * vrhs[i][j];
        }
    }
    set_boundary_condition_V(u, v, wing_x, wing_y);  // 境界条件を適用
}


/* csv入力 */
std::vector<int> read_csv(std::string fileName) {
    std::vector<int> data;
    std::ifstream ifs(fileName);
    std::string str;

    if (!ifs) {
		std::cout << "! couldn't open .csv file !\n";
		exit(0);
    }

    while (getline(ifs, str)) { //一行ずつ抽出
        std::string token;
        std::istringstream iss(str);
        while (getline(iss, token, ',')) {  //一行内をカンマで分解
            data.push_back(std::stof(token.c_str()));   //文字列データなので数値に変換
        }
    }
    return data;
}

/* csv出力 */
int write_csv_puv(double puv[mx + 2][my + 2], const char *file_name) {
    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = my; i >= 1; i--) {
        for (int j = 1; j <= mx; j++) {
            fprintf(fp, "%f", puv[j][i]);
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
int write_csv_cdclcp(double cd, double cl, double cp1, double cp2, const char *file_name) {
    FILE *fp;
    fp = fopen(file_name, "a");
    fprintf(fp, "%f,%f,%f,%f\n", cd, cl, cp1, cp2);
    fclose(fp);
    return 0;
}



int main() {
    const std::string alpha("20"); // 適当な迎角を設定

    // 翼の座標をcsvから得る
    std::vector<int> data1 = read_csv("./airfoil/wing_data/joukowsky_" + alpha + ".csv");
    std::vector<int> wing_x, wing_y;
    for (int i=0; i<(int)data1.size()/2; i++) {
        wing_x.push_back(data1[i]);
        wing_y.push_back(data1[i+data1.size()/2]);
    }
    std::vector<int> data2 = read_csv("./airfoil/wing_data/joukowsky_wall_" + alpha + ".csv");
    std::vector<int> wing_wall_x, wing_wall_y;
    for (int i=0; i<(int)data2.size()/2; i++) {
        wing_wall_x.push_back(data2[i]);
        wing_wall_y.push_back(data2[i+data2.size()/2]);
    }

    // 初期条件をセット
    set_init_condition(u, v, p);
    set_boundary_condition_p(p, wing_wall_x, wing_wall_y, wing_x, wing_y);
    set_boundary_condition_V(u, v, wing_x, wing_y);

    // write_csv_cdclcp()の準備
    std::string s_cdclcp = "./airfoil/output_csv/cdclcp_" + alpha + ".csv";
    const char *cs_cdclcp = s_cdclcp.data();


    std::chrono::system_clock::time_point start_time, end_time;  // 時間計測用変数を用意
    start_time = std::chrono::system_clock::now();  // 計測始まり

    // メインルーチン
    for (int n = 1; n <= nlast; n++) {
		solve_poisson_eq(u, v, p, dx, dy, dt, wing_wall_x, wing_wall_y, wing_x, wing_y);
		solve_velocity_eq(u, v, p, dx, dy, dt, wing_wall_x, wing_wall_y, wing_x, wing_y);

		for (int i = 1; i <= mx; i++) {
			for (int j = 1; j <= my; j++) {
				cp[i][j] = 2 * p[i][j];
			}
		}

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
		double cp1 = 2 * p[3 * i_2 - 2 * i_1][j_1];
		double cp2 = 2 * p[3 * i_2 - 2 * i_1][j_2];

		end_time = std::chrono::system_clock::now();  // 計測終わり
		double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time -start_time).count() / 1000.0;	 // かかった時間[s]
		std::cout << "step:" << n << "/" << nlast << " elapsed_time:" << elapsed_time << "[s]" << std::endl;
        
        // 結果の出力
        if (n % 50 == 0) {
            std::string s_p = "./airfoil/output_csv/p_" + alpha + "_" + std::to_string(n) + ".csv";
            const char *cs_p = s_p.data();
            write_csv_puv(cp, cs_p);
            std::string s_u = "./airfoil/output_csv/u_" + alpha + "_" + std::to_string(n) + ".csv";
            const char *cs_u = s_u.data();
            write_csv_puv(u, cs_u);
            std::string s_v = "./airfoil/output_csv/v_" + alpha + "_" + std::to_string(n) + ".csv";
            const char *cs_v = s_v.data();
            write_csv_puv(v, cs_v);
        }
		write_csv_cdclcp(cd, cl, cp1, cp2, cs_cdclcp);
    }

    return 0;
}

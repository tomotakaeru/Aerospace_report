#include<stdio.h>
#include<cmath>
#include<vector>
#include<fstream>
#include<sstream>

//条件の設定
const double rho = 1.225;   //密度[kg/m^3]
const double T0 = 293.0;    //温度[K]
const double mu = 1.458e-06 * std::pow(T0, 1.5) / (T0+110.4);    //動粘性係数[kg/m*s]
const double width = 0.5;    //x範囲[m]
const double U0 = 20.0; //x方向初期速度[m/s]
const double Re = rho * U0 * width / mu; //レイノルズ数
const double height = 2 * 5.2 * width / std::sqrt(Re);   //y範囲[m],Blausiusの解の2倍
const double BL_threshold = 0.995;  //境界層厚さの基準

//計算格子
const int NX = 50000;
const int NY = 250;
const double dx = width / NX;
const double dy = height / NY;
double u[NY+1][NX+1];
double v[NY+1][NX+1];



int main()
{
    //初期条件
    for (int j=0; j<NY+1; j++){
        u[j][0] = U0;
        v[j][0] = 0.0;
    }

    //i列のデータを一時保存する配列
    double u_pre[NY+1];
    double v_pre[NY+1];

    for (int i=0; i<NX; i++){
        //i列のデータを一時保存
        for (int j=0; j<NY+1; j++){
            u_pre[j] = u[j][i];
            v_pre[j] = v[j][i];
        }
        //i+1列のuを境界除いて求める(運動量保存則)
        for (int j=1; j<NY; j++){
            double dudy = (mu/rho * (u_pre[j+1]-2*u_pre[j]+u_pre[j-1])/(dy*dy) - v_pre[j] * (u_pre[j+1]-u_pre[j-1])/(2*dy));
            u[j][i+1] = u_pre[j] + dudy*dx / u_pre[j];
        }
        //i+1列のu,vの境界条件
        u[0][i+1] = 0.0;
        u[NY][i+1] = U0;
        if (i*dx>=0.2 && i*dx<=0.3){
            v[0][i+1] = -0.04;
        }else{
            v[0][i+1] = 0.0;
        }
        //i+1列のvを境界除いて求める(質量保存則)
        for (int j=1; j<NY+1; j++){
            double dudx = (u[j][i+1]-u_pre[j])/dx + (u[j+1][i+1]-u_pre[j+1])/dx;
            v[j][i+1] = v[j-1][i+1] - dudx*dy/2;
        }
    }

    //u,vをcsv出力
    int split_num = 10;
    for (int i = 0; i < split_num+1; i++) {
        int output_cols = NX / split_num * i;
        std::stringstream ss1;
        ss1 << "./div_minus/uv_x" << output_cols*dx << ".csv";
        std::string filename1 = ss1.str();
        std::ofstream ofs1(filename1);
        ofs1 << "y,u,v" << std::endl;
        for (int j=0; j<NY+1; j++){
            ofs1 << j*dy*1000 << "," << u[j][output_cols] << "," << v[j][output_cols] << std::endl;
        }
	}


    //これ以降のcsv出力準備
    std::stringstream ss2;
    ss2 << "./div_minus/thickness.csv";
    std::string filename2 = ss2.str();
    std::ofstream ofs2(filename2);
    for (int i=0; i<NX+1; i++){
        if (i==NX){
            ofs2 << dx*i << std::endl;
        }else{
            ofs2 << dx*i << ",";
        }
    }

    //境界層厚さを求めて出力
    std::vector<double> delta;
    for (int i=0; i<NX+1; i++){
        int j=0;
        while (u[j][i] < U0*BL_threshold){
            j++;
        }
        delta.push_back(j*dy*1000);
        if (i==NX){
            ofs2 << delta[i] << std::endl;
        }else{
            ofs2 << delta[i] << ",";
        }
    }

    //排除厚さを求めて出力
    std::vector<double> delta_star;
    for (int i=0; i<NX+1; i++){
        double delta_star_pre=0;
        int j=0;
        while (std::abs( 1 - u[j][i]/U0 ) > 0.0001){
            delta_star_pre += (1 - u[j][i]/U0) * dy;
            j++;
        }
        delta_star.push_back(delta_star_pre*1000);
        if (i==NX){
            ofs2 << delta_star[i] << std::endl;
        }else{
            ofs2 << delta_star[i] << ",";
        }
    }

    //運動量厚さを求めて出力
    std::vector<double> theta;
    for (int i=0; i<NX+1; i++){
        double theta_pre=0;
        int j=0;
        while (std::abs( 1 - u[j][i]/U0 ) > 0.0001){
            theta_pre += u[j][i]/U0 * (1 - u[j][i]/U0) * dy;
            j++;
        }
        theta.push_back(theta_pre*1000);
        if (i==NX){
            ofs2 << theta[i] << std::endl;
        }else{
            ofs2 << theta[i] << ",";
        }
    }

    //エネルギー厚さを求めて出力
    std::vector<double> Theta;
    for (int i=0; i<NX+1; i++){
        double Theta_pre=0;
        int j=0;
        while (std::abs( 1 - std::pow(u[j][i]/U0, 2) ) > 0.0001){
            Theta_pre += u[j][i]/U0 * (1 - std::pow(u[j][i]/U0, 2)) * dy;
            j++;
        }
        Theta.push_back(Theta_pre*1000);
        if (i==NX){
            ofs2 << Theta[i] << std::endl;
        }else{
            ofs2 << Theta[i] << ",";
        }
    }

    //壁面摩擦係数を求めて出力
    std::vector<double> Cf;
    for (int i=0; i<NX+1; i++){
        Cf.push_back(2/(rho*U0*U0) * (mu * (u[1][i]-u[0][i])/dy) * 1000);
        if (i==NX){
            ofs2 << Cf[i] << std::endl;
        }else{
            ofs2 << Cf[i] << ",";
        }
    }



    //print
    printf("流速：\t\t\t%.3f[m/s]\n", U0);
    printf("レイノルズ数：\t\t%5.3e\n", Re);
    printf("x方向範囲：\t\t%.3f[m]\n", width);
    printf("y方向範囲：\t\t%.3f[mm]\n", height*1000);
    printf("x方向分割数：\t\t%d\n", NX);
    printf("y方向分割数：\t\t%d\n", NY);
    printf("境界層厚さ（厳密解）：\t%f[mm]\n", height/2*1000);
    printf("境界層厚さ：\t\t%f[mm]\n", delta[NX]);
    printf("排除厚さ：\t\t%f[mm]\n", delta_star[NX]);
    printf("運動量厚さ：\t\t%f[mm]\n", theta[NX]);
    printf("エネルギー厚さ：\t%f[mm]\n", Theta[NX]);
    printf("壁面摩擦係数：\t\t%f\n", Cf[NX]);

    return 0;
}

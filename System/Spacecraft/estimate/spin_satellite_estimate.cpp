#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <random>

using namespace Eigen;


Vector3d solve_omega_eq(Vector3d omega, Matrix3d I, Vector3d M){
    PartialPivLU<Matrix3d> lu(I);
    return lu.solve(M - omega.cross(I * omega));
}
Vector3d RungeKutta_omega(Vector3d omega, Matrix3d I, Vector3d M, double dt){
    Vector3d k0 = solve_omega_eq(omega, I, M) * dt;
    Vector3d k1 = solve_omega_eq(omega + k0 / 2, I, M) * dt;
    Vector3d k2 = solve_omega_eq(omega + k1 / 2, I, M) * dt;
    Vector3d k3 = solve_omega_eq(omega + k2, I, M) * dt;
    return (k0 + k1 * 2 + k2 * 2 + k3) / 6;
}
Vector4d solve_q_eq(Vector4d q, Vector3d omega){
    Matrix<double, 4, 3> Q;
    Q << -q(1), -q(2), -q(3),
          q(0), -q(3),  q(2),
          q(3),  q(0), -q(1),
         -q(2),  q(1),  q(0);
    return Q * omega / 2;
}
Quaterniond RungeKutta_q(Quaterniond q, Vector3d omega, Vector3d omega_half_next, Vector3d omega_next, double dt){
    Vector4d q_vec(q.w(), q.x(), q.y(), q.z());
    Vector4d k0 = solve_q_eq(q_vec, omega) * dt;
    Vector4d k1 = solve_q_eq(q_vec + k0 / 2, omega_half_next) * dt;
    Vector4d k2 = solve_q_eq(q_vec + k1 / 2, omega_half_next) * dt;
    Vector4d k3 = solve_q_eq(q_vec + k2, omega_next) * dt;
    Vector4d tmp = (k0 + k1 * 2 + k2 * 2 + k3) / 6;
    return Quaterniond (tmp(0), tmp(1), tmp(2), tmp(3));
}

Matrix<double, 7, 7> calc_A(Quaterniond q, Vector3d omega, Matrix3d I){
    Matrix<double, 7, 7> result;
    result <<        0, -0.5 * omega(0), -0.5 * omega(1), -0.5 * omega(2), -0.5 * q.x(), -0.5 * q.y(), -0.5 * q.z(),
        0.5 * omega(0),               0,  0.5 * omega(2), -0.5 * omega(1),  0.5 * q.w(), -0.5 * q.z(),  0.5 * q.y(),
        0.5 * omega(1), -0.5 * omega(2),               0,  0.5 * omega(0),  0.5 * q.z(),  0.5 * q.w(), -0.5 * q.x(),
        0.5 * omega(2),  0.5 * omega(1), -0.5 * omega(0),               0, -0.5 * q.y(),  0.5 * q.x(),  0.5 * q.w(),
                     0,               0,               0,               0, 0, (I(1,1) - I(2,2)) / I(0,0) * omega(2), (I(1,1) - I(2,2)) / I(0,0) * omega(1),
                     0,               0,               0,               0, (I(2,2) - I(0,0)) / I(1,1) * omega(2), 0, (I(2,2) - I(0,0)) / I(1,1) * omega(0),
                     0,               0,               0,               0, (I(0,0) - I(1,1)) / I(2,2) * omega(1), (I(0,0) - I(1,1)) / I(2,2) * omega(0), 0;
    return result;
}
Matrix<double, 3, 7> calc_H(Quaterniond q, int DCM_index){
    Matrix<double, 3, 7> result;
    switch(DCM_index){
    case 0:
        result << 2 * q.w(), 2 * q.x(), -2 * q.y(), -2 * q.z(), 0, 0, 0,
                  2 * q.z(), 2 * q.y(),  2 * q.x(),  2 * q.w(), 0, 0, 0,
                 -2 * q.y(), 2 * q.z(), -2 * q.w(),  2 * q.x(), 0, 0, 0;
        break;
    case 1:
        result << -2 * q.z(),  2 * q.y(), 2 * q.x(), -2 * q.w(), 0, 0, 0,
                   2 * q.w(), -2 * q.x(), 2 * q.y(), -2 * q.z(), 0, 0, 0,
                   2 * q.x(),  2 * q.w(), 2 * q.z(),  2 * q.y(), 0, 0, 0;
        break;
    case 2:
        result << 2 * q.y(),  2 * q.z(),  2 * q.w(), 2 * q.x(), 0, 0, 0,
                 -2 * q.x(), -2 * q.w(),  2 * q.z(), 2 * q.y(), 0, 0, 0,
                  2 * q.w(), -2 * q.x(), -2 * q.y(), 2 * q.z(), 0, 0, 0;
        break;
    }
    return result;
}


int main()
{
    /* 計算条件 */
    const double dt = 0.01;    //[s]
    const double dt_obs = 1;
    const double end_time = 50;
    const int N = (int)(end_time / dt);
    const int N_obs = (int)(dt_obs / dt);
    Matrix3d I;  //[kgm^2]
    I << 1.9, 0.0, 0.0,
         0.0, 1.6, 0.0,
         0.0, 0.0, 2.0;
    Vector3d M_O(0.0, 0.0, 0.0);  //[Nm]
    Vector3d M_C(0.0, 0.0, 0.0);
    const double omega_s = 17 * 2 * M_PI / 60;  //17[rpm]to[rad/s]
    const Vector3d omega_0(0.1, omega_s + 0.1, 0.0);  //[rad/s]
    const Quaterniond q_0(1.0, 0.0, 0.0, 0.0);
    const Vector3d domega_0(0.1, -0.1, 0.1);
    const Quaterniond dq_0(0.1, -0.1, 0.1, -0.1);
    const double sigma_w = 0.01;
    const double sigma_v = 0.01;

    /* 初期化 */
    Vector3d omega = omega_0;
    Quaterniond q = q_0;
    Vector3d omega_hat = omega_0 + domega_0;
    Quaterniond q_hat = q_0;
    q_hat.coeffs() += dq_0.coeffs();
    Vector3d omega_half_next;
    Vector3d omega_hat_half_next;
    Vector3d omega_next;
    Vector3d omega_hat_next;
    Vector3d DCM_row;
    Vector3d DCM_row_hat;
    Vector3d DCM_err;
    Vector3d::Index min_index;
    int DCM_index;
    Vector3d v;
    Matrix<double, 7, 7> A;
    Matrix<double, 7, 7> A_k;
    Matrix<double, 7, 3> B = MatrixXd::Zero(7, 3);
    B(4,0) = 1.0 / I(0,0);
    B(5,1) = 1.0 / I(1,1);
    B(6,2) = 1.0 / I(2,2);
    Matrix<double, 7, 3> B_k;
    Matrix<double, 3, 7> H;
    Matrix<double, 7, 3> K;
    Matrix<double, 7, 7> P = sigma_w * MatrixXd::Identity(7, 7);
    Matrix<double, 7, 7> M = sigma_w * MatrixXd::Identity(7, 7);
    Matrix<double, 3, 3> Q = sigma_w*sigma_w * Matrix3d::Identity(3, 3);
    Matrix<double, 3, 3> R = sigma_v*sigma_v * Matrix3d::Identity(3, 3);
    VectorXd x_hat(7);

    /* DCM選択用＆ノイズ用の乱数生成の準備 */
    // std::random_device seed; //非決定的乱数生成の場合
    // std::mt19937 mt(seed());
    const int seed = 9; //seed固定の場合
    std::mt19937 mt(seed);
    std::uniform_int_distribution<> select_DCM(0, 2);
    std::normal_distribution<> noise_w(0.0, sigma_w*sigma_w); 
    std::normal_distribution<> noise_v(0.0, sigma_v*sigma_v); 

    /* csvへの結果出力の準備 */
    std::ofstream ofs;
    ofs.open("./2/output_omega_q.csv", std::ios::trunc);
    ofs << "omega_x, omega_y, omega_z, q_0, q_1, q_2, q_3, omega_x_hat, omega_y_hat, omega_z_hat, q_0_hat, q_1_hat, q_2_hat, q_3_hat, P11, P22, P33, P44, P55, P66, P77" << std::endl;
    ofs << omega(0) << "," << omega(1) << "," << omega(2) << "," << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << "," << omega_hat(0) << "," << omega_hat(1) << "," << omega_hat(2) << "," << q_hat.w() << "," << q_hat.x() << "," << q_hat.y() << "," << q_hat.z() << "," << P(0,0) << "," << P(1,1) << "," << P(2,2) << "," << P(3,3) << "," << P(4,4) << "," << P(5,5) << "," << P(6,6) << std::endl;

    /* メインルーチン */
    for(int i=1; i<N+1; i++){
        /* 状態量真値更新 */
        M_O << noise_w(mt), noise_w(mt), noise_w(mt);
        omega_half_next = omega + RungeKutta_omega(omega, I, M_O + M_C, dt/2);
        M_O << noise_w(mt), noise_w(mt), noise_w(mt);
        omega_next = omega_half_next + RungeKutta_omega(omega_half_next, I, M_O + M_C, dt/2);
        q.coeffs() += RungeKutta_q(q, omega, omega_half_next, omega_next, dt).coeffs();
        q.normalize();
        omega = omega_next;
        /* 状態量推定値更新 */
        M_O << 0, 0, 0;
        omega_hat_half_next = omega_hat + RungeKutta_omega(omega_hat, I, M_O + M_C, dt/2);
        omega_hat_next = omega_hat_half_next + RungeKutta_omega(omega_hat_half_next, I, M_O + M_C, dt/2);
        q_hat.coeffs() += RungeKutta_q(q_hat, omega_hat, omega_hat_half_next, omega_hat_next, dt).coeffs();
        q_hat.normalize();
        omega_hat = omega_hat_next;
        /* 状態量推定値共分散更新 */
        A = calc_A(q_hat, omega_hat, I);
        A_k = A * dt + MatrixXd::Identity(7, 7);
        B_k = B * dt;
        // A_k = (A * dt).array().exp().matrix();
        // B_k = A.inverse() * (A_k - MatrixXd::Identity(7, 7)) * B;
        M = A_k * P * A_k.transpose() + B_k * Q * B_k.transpose();
        P = M;

        if(i % N_obs == 0){ //観測を行うタイミング
            // DCM_index = 1;
            // DCM_index = select_DCM(mt);
            // DCM_row = q.toRotationMatrix().block(0, DCM_index, 3, 1);
            // v << noise_v(mt), noise_v(mt), noise_v(mt);
            // DCM_row += v;
            // DCM_row_hat = q_hat.toRotationMatrix().block(0, DCM_index, 3, 1);
            /* 観測量真値生成 */
            DCM_row = q.toRotationMatrix().block(0, select_DCM(mt), 3, 1);
            v << noise_v(mt), noise_v(mt), noise_v(mt);
            DCM_row += v;
            /* 観測量推定値を選択 */
            for(int j = 0;j < 3;j++){
                DCM_err(j) = (DCM_row - q_hat.toRotationMatrix().block(0, j, 3, 1)).norm();
            }
            DCM_err.minCoeff(&min_index);
            DCM_index = min_index;
            DCM_row_hat = q_hat.toRotationMatrix().block(0, DCM_index, 3, 1);
            /* 状態量推定値とその共分散を更新 */
            H = calc_H(q_hat, DCM_index); 
            P = M - M * H.transpose() * (H * M * H.transpose() + R).inverse() * H * M;
            K = P * H.transpose() * R.inverse();
            x_hat = K * (DCM_row - DCM_row_hat);
            Quaterniond dq(x_hat(0), x_hat(1), x_hat(2), x_hat(3));  //q.coeffs()=(q1,q2,q3,q0)の順のため
            q_hat.coeffs() += dq.coeffs();
            q_hat.normalize();
            omega_hat += x_hat.tail(3);
            std::cout << "q_hat_after:\n" << q_hat.coeffs() << std::endl;
        }
        ofs << omega(0) << "," << omega(1) << "," << omega(2) << "," << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << "," << omega_hat(0) << "," << omega_hat(1) << "," << omega_hat(2) << "," << q_hat.w() << "," << q_hat.x() << "," << q_hat.y() << "," << q_hat.z() << "," << P(0,0) << "," << P(1,1) << "," << P(2,2) << "," << P(3,3) << "," << P(4,4) << "," << P(5,5) << "," << P(6,6) << std::endl;
        std::cout << "output" << i << std::endl;
    }
    ofs.close();
}

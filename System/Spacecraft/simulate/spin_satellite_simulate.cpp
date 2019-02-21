#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

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


int main()
{
    /* 計算条件 */
    const double dt = 0.01;    //[s]
    const double end_time = 100;
    const int N = (int)(end_time / dt);
    Matrix3d I;  //[kgm^2]
    I << 1.9, 0.0, 0.0,
         0.0, 1.6, 0.0,
         0.0, 0.0, 2.0;
    Vector3d M_O(0.0, 0.0, 0.0);  //[Nm]
    Vector3d M_C(0.0, 0.0, 0.0);
    const double omega_s = 17 * 2 * M_PI / 60;  //17[rpm]to[rad/s]
    const Vector3d omega_0(0.1, omega_s + 0.1, 0.0);  //[rad/s]
    const Quaterniond q_0(1.0, 0.0, 0.0, 0.0);

    /* 初期化 */
    Vector3d omega = omega_0;
    Quaterniond q = q_0;
    Matrix3d DCM = q.toRotationMatrix();
    Vector3d omega_half_next;
    Vector3d omega_next;

    /* csvへの結果出力の準備 */
    std::ofstream ofs1;
    ofs1.open("./1/output_omega_q.csv", std::ios::trunc);
    ofs1 << "omega_x, omega_y, omega_z, q_0, q_1, q_2, q_3" << std::endl;
    ofs1 << omega(0) << "," << omega(1) << "," << omega(2) << "," << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << std::endl;
    std::ofstream ofs2;
    ofs2.open("./1/output_DCM.csv", std::ios::trunc);
    ofs2 << "11, 21, 31, 12, 22, 32, 13, 23, 33" << std::endl;
    ofs2 << DCM(0,0) << "," << DCM(0,1) << "," << DCM(0,2) << "," << DCM(1,0) << "," << DCM(1,1) << "," << DCM(1,2) << "," << DCM(2,0) << "," << DCM(2,1) << "," << DCM(2,2) << std::endl;

    /* メインルーチン */
    for(int i=0; i<N; i++){
        omega_half_next = omega + RungeKutta_omega(omega, I, M_O + M_C, dt/2);
        omega_next = omega_half_next + RungeKutta_omega(omega_half_next, I, M_O + M_C, dt/2);
        q.coeffs() += RungeKutta_q(q, omega, omega_half_next, omega_next, dt).coeffs();
        q.normalize();
        DCM = q.toRotationMatrix();
        omega = omega_next;

        ofs1 << omega(0) << "," << omega(1) << "," << omega(2) << "," << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << std::endl;
        ofs2 << DCM(0,0) << "," << DCM(0,1) << "," << DCM(0,2) << "," << DCM(1,0) << "," << DCM(1,1) << "," << DCM(1,2) << "," << DCM(2,0) << "," << DCM(2,1) << "," << DCM(2,2) << std::endl;
        std::cout << "output" << i+1 << std::endl;
    }
    ofs1.close();
    ofs2.close();
}

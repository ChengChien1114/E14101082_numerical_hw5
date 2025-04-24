#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// --- 問題 1: y' = 1 + (y/t) + (y/t)^2, y(1) = 0 ---
double f1(double t, double y) {
    if (t == 0) { 
        return 1.0; 
    }
    double term = y / t;
    return 1.0 + term + term * term;
}

// Exact solution for Problem 1: y(t) = t * tan(ln(t))
double y_exact(double t) {
    if (t <= 0) return 0;
    if (t == 1.0) return 0.0; 
    return t * tan(log(t)); 
}

// --- 問題 1b: Taylor Order 2 ---
double df_dt(double t, double y) {
    if (t == 0) return 0;
    return -y / (t * t) - 2.0 * y * y / (t * t * t);
}

double df_dy(double t, double y) {
    if (t == 0) return 0;
    return 1.0 / t + 2.0 * y / (t * t);
}

double f1_prime(double t, double y) {
    if (t == 0) return 0;
    return df_dt(t, y) + f1(t, y) * df_dy(t, y);
}


int main() {
    std::cout << std::fixed << std::setprecision(8);

    // --- 問題 1 求解 ---
    std::cout << "--- Problem 1: y' = 1+(y/t)+(y/t)^2, y(1)=0, 1 <= t <= 2 ---" << std::endl;
    double t0_1 = 1.0, y0_1 = 0.0;
    double t_end_1 = 2.0;
    double h_1 = 0.1;
    int n_steps_1 = static_cast<int>(std::round((t_end_1 - t0_1) / h_1));

    // a) Euler's Method
    std::cout << "\n--- Part (a): Euler's Method (h = " << h_1 << ") ---" << std::endl;
    std::cout << "t_i       w_i (Euler)    y_i (Exact)    |Error|" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    double t_euler = t0_1;
    double w_euler = y0_1;
    std::cout << t_euler << "  " << w_euler << "  " << y_exact(t_euler) << "  " << std::abs(w_euler - y_exact(t_euler)) << std::endl;
    for (int i = 0; i < n_steps_1; ++i) {
        double current_h = std::min(h_1, t_end_1 - t_euler);
        if (current_h <= 1e-9) break; 
        w_euler = w_euler + current_h * f1(t_euler, w_euler);
        t_euler = t_euler + current_h; 
        double t_print = std::round(t_euler * 100.0) / 100.0;

        std::cout << t_print << "  " << w_euler << "  " << y_exact(t_euler) << "  " << std::abs(w_euler - y_exact(t_euler)) << std::endl;
    }

    // b) Taylor's Method Order 2
    std::cout << "\n--- Part (b): Taylor's Method Order 2 (h = " << h_1 << ") ---" << std::endl;
    std::cout << "t_i       w_i (Taylor2)  y_i (Exact)    |Error|" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    double t_taylor2 = t0_1;
    double w_taylor2 = y0_1;
    std::cout << t_taylor2 << "  " << w_taylor2 << "  " << y_exact(t_taylor2) << "  " << std::abs(w_taylor2 - y_exact(t_taylor2)) << std::endl;
    for (int i = 0; i < n_steps_1; ++i) {
        double current_h = std::min(h_1, t_end_1 - t_taylor2);
         if (current_h <= 1e-9) break;
        double T2 = f1(t_taylor2, w_taylor2) + (current_h / 2.0) * f1_prime(t_taylor2, w_taylor2);
        w_taylor2 = w_taylor2 + current_h * T2;
        t_taylor2 = t_taylor2 + current_h; 
        double t_print = std::round(t_taylor2 * 100.0) / 100.0;

        std::cout << t_print << "  " << w_taylor2 << "  " << y_exact(t_taylor2) << "  " << std::abs(w_taylor2 - y_exact(t_taylor2)) << std::endl;
    }

    return 0;
}
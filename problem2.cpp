#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> 

// --- 問題 2: System u1', u2' ---
// f1(t, u1, u2) for Problem 2: u1' = 9*u1 + 24*u2 + 5*cos(t) - (1/3)*sin(t)
double f2_1(double t, double u1, double u2) {
    return 9.0 * u1 + 24.0 * u2 + 5.0 * cos(t) - (1.0 / 3.0) * sin(t);
}

// f2(t, u1, u2) for Problem 2: u2' = -24*u1 - 51*u2 - 9*cos(t) + (1/3)*sin(t)
double f2_2(double t, double u1, double u2) {
    return -24.0 * u1 - 51.0 * u2 - 9.0 * cos(t) + (1.0 / 3.0) * sin(t);
}

// Exact solution u1(t) for Problem 2
double u1_exact(double t) {
    return 2.0 * exp(-3.0 * t) - exp(-39.0 * t) + (1.0 / 3.0) * cos(t);
}

// Exact solution u2(t) for Problem 2
double u2_exact(double t) {
    return -exp(-3.0 * t) + 2.0 * exp(-39.0 * t) - (1.0 / 3.0) * cos(t);
}


int main() {
    std::cout << std::fixed << std::setprecision(8);
    // --- 問題 2 求解 ---
    std::cout << "\n--- Problem 2: System u1', u2', u1(0)=4/3, u2(0)=2/3 ---" << std::endl;
    double t0_2 = 0.0;
    double u1_0 = 4.0 / 3.0;
    double u2_0 = 2.0 / 3.0;
    double t_end_2 = 1.0;

    for (double h_2 : {0.1, 0.05}) { 
         std::cout << "\n--- Runge-Kutta Order 4 (h = " << h_2 << ", t_end = " << t_end_2 << ") ---" << std::endl;
         std::cout << "t_j       w1_j (RK4)   u1_j (Exact) |Error u1|   w2_j (RK4)   u2_j (Exact) |Error u2|" << std::endl;
         std::cout << "------------------------------------------------------------------------------------------" << std::endl;
         double t_rk4 = t0_2;
         double w1_rk4 = u1_0;
         double w2_rk4 = u2_0;
         int n_steps_2 = static_cast<int>(std::round((t_end_2 - t0_2) / h_2));


         std::cout << t_rk4 << "  "
                   << w1_rk4 << "  " << u1_exact(t_rk4) << "  " << std::abs(w1_rk4 - u1_exact(t_rk4)) << "   "
                   << w2_rk4 << "  " << u2_exact(t_rk4) << "  " << std::abs(w2_rk4 - u2_exact(t_rk4)) << std::endl;

         for (int j = 0; j < n_steps_2; ++j) {
             double current_h = std::min(h_2, t_end_2 - t_rk4);
             if (current_h <= 1e-9) break; 
             // Calculate k1
             double k1_1 = current_h * f2_1(t_rk4, w1_rk4, w2_rk4);
             double k1_2 = current_h * f2_2(t_rk4, w1_rk4, w2_rk4);
             // Calculate k2
             double k2_1 = current_h * f2_1(t_rk4 + current_h / 2.0, w1_rk4 + k1_1 / 2.0, w2_rk4 + k1_2 / 2.0);
             double k2_2 = current_h * f2_2(t_rk4 + current_h / 2.0, w1_rk4 + k1_1 / 2.0, w2_rk4 + k1_2 / 2.0);
             // Calculate k3
             double k3_1 = current_h * f2_1(t_rk4 + current_h / 2.0, w1_rk4 + k2_1 / 2.0, w2_rk4 + k2_2 / 2.0);
             double k3_2 = current_h * f2_2(t_rk4 + current_h / 2.0, w1_rk4 + k2_1 / 2.0, w2_rk4 + k2_2 / 2.0);
             // Calculate k4
             double k4_1 = current_h * f2_1(t_rk4 + current_h, w1_rk4 + k3_1, w2_rk4 + k3_2);
             double k4_2 = current_h * f2_2(t_rk4 + current_h, w1_rk4 + k3_1, w2_rk4 + k3_2);
             // Update w
             w1_rk4 = w1_rk4 + (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1) / 6.0;
             w2_rk4 = w2_rk4 + (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2) / 6.0;
             t_rk4 = t_rk4 + current_h; 
             double t_print = std::round(t_rk4 * 1000.0) / 1000.0; 
             std::cout << t_print << "  "
                       << w1_rk4 << "  " << u1_exact(t_rk4) << "  " << std::abs(w1_rk4 - u1_exact(t_rk4)) << "   "
                       << w2_rk4 << "  " << u2_exact(t_rk4) << "  " << std::abs(w2_rk4 - u2_exact(t_rk4)) << std::endl;
         }
    }

    return 0;
}
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


unsigned long long n = 2; // system size
const double a = 3;       // variant parameter
double epsilon;           // error
unsigned long long iteration_cnt = 0;

std::ofstream answer_file;


unsigned long long find_max(std::vector<std::vector<double>> a, unsigned long long k) {
    double max = std::fabs(a[k][k]);
    unsigned long long max_id = k;
    for (unsigned long long i = k+1; i < n; ++i) { // index of a row
        if (std::fabs(a[i][k]) > max) {
            max = std::fabs(a[i][k]);
            max_id = i;
        }
    }

    return max_id; 
}

void swap_rows(std::vector<std::vector<double>> &a, unsigned long long k, unsigned long long max_id) {
    for (unsigned long long j = 0; j < n; ++j) // index of a column
        std::swap(a[k][j], a[max_id][j]);
}

void swap_columns(std::vector<std::vector<double>> &a, unsigned long long k, unsigned long long max_id) {
    for (unsigned long long i = 0; i < n; ++i) // index of a row
        std::swap(a[i][k], a[i][max_id]);
}

void make_identity_matrix(std::vector<std::vector<double>> &I) {
    for (unsigned long long i = 0; i < n; ++i)
        I[i][i] = 1; 
}

void LU_decomposition(std::vector<std::vector<double>> a, std::vector<double> &B, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U) {
    double mu;

    unsigned long long max_id;
    make_identity_matrix(L);
    for (unsigned long long k = 0; k < (n - 1); ++k) { // step
        max_id = find_max(a, k);
        if (max_id != k) {
            swap_rows(a, k, max_id);
            std::swap(B[k], B[max_id]);
            swap_rows(L, k, max_id);
            swap_columns(L, k, max_id);
        }

        for (unsigned long long i = (k + 1); i < n; ++i) { // index of a row
            mu = a[i][k] / a[k][k]; // coefficient to make k-th element of current row equal to zero
            L[i][k] = mu;
            for (unsigned long long j = k; j < n; ++j) { // index of a column
                a[i][j] -= mu * a[k][j]; // subtract k-th row multiplied by mu from all rows below it
            }
        }        
    }
    L[n-1][n-1] = 1;
    U = a;
}

// solves Lz = b 
std::vector<double> find_z(std::vector<std::vector<double>> l, std::vector<double> &b) {
    std::vector<double> z(n);
    z[0] = b[0];
    for (unsigned long long i = 1; i < n; ++i) { // index of a row
        for (unsigned long long j = 0; j < i; ++j) { // index of a column
            b[i] -= l[i][j] * z[j];
        }
        z[i] = b[i];
    } 

    return z;
}

// solves Ux = z 
std::vector<double> find_x(std::vector<std::vector<double>> u, std::vector<double> z) {
    std::vector<double> x(n);
    for (unsigned long long i = n; i-- > 0;) { // index of a row
        for (unsigned long long j = n; j-- > (i + 1);) { // index of a column
            z[i] -= u[i][j] * x[j];
        }
        x[i] = z[i] / u[i][i];
    } 

    return x;
}


double phi1(double x2) {
    return std::cos(x2) + a;
}

double phi1_d(double x2) {
    return -std::sin(x2);
}

double phi2(double x1) {
    return std::sin(x1) + a;
}

double phi2_d(double x1) {
    return std::cos(x1);
}

double J_phi(double x1, double x2) {
    return -(phi1_d(x2) * phi2_d(x1));
}


std::pair<double, double> fixed_point_iteration(double l1, double r1, double l2, double r2) {
    iteration_cnt = 0;
    double x1_cur;
    double x2_cur;
    double x1_prev = r1;
    double x2_prev = r2;

    double q = std::max(std::max(std::abs(J_phi(l1, r1)), 
                                 std::abs(J_phi(l1, r2))), 
                        std::max(std::abs(J_phi(l2, r1)), 
                                 std::abs(J_phi(l2, r2))));
    double epsilon_k;
    do {
        x1_cur = phi1(x2_prev);
        x2_cur = phi2(x1_prev);
        iteration_cnt++;
        epsilon_k = (q / (1 - q)) * (std::abs(x1_cur - x1_prev) + std::abs(x2_cur - x2_prev));
        x1_prev = x1_cur;
        x2_prev = x2_cur;
    } while (epsilon_k > epsilon);

    return std::make_pair(x1_cur, x2_cur);
}


double f1(double x1, double x2) {
    return x1 - std::cos(x2) - a;
}

double f2(double x1, double x2) {
    return x2 - std::sin(x1) - a;
}

std::vector<std::vector<double>> J_f(double x1, double x2) {
    std::vector<std::vector<double>> J_f(n);
    for (unsigned long long i = 0; i < n; ++i)
        J_f[i].resize(n);

    J_f[0][0] = 1.0;            // df1/dx1
    J_f[0][1] = std::sin(x2);   // df1/dx2
    J_f[1][0] = -std::cos(x1);  // df2/dx1
    J_f[1][1] = 1.0;            // df2/dx2
    return J_f;
}

std::pair<double, double> newtons_method(double x1_0, double x2_0) {
    iteration_cnt = 0;
    std::vector<double> x_cur(n);
    std::vector<double> x_prev = {x1_0, x2_0};

    double x1 , x2;
    std::vector<double> minus_f;
    std::vector<double> delta_x(n);
    double epsilon_k;
    do {
        x1 = x_prev[0];
        x2 = x_prev[1];
        minus_f = {-f1(x1, x2), -f2(x1, x2)};

        std::vector<std::vector<double>> L(n);
        std::vector<std::vector<double>> U(n);
        for (unsigned long long i = 0; i < n; ++i) {
            L[i].resize(n);
            U[i].resize(n);
        }
        LU_decomposition(J_f(x1, x2), minus_f, L, U); 
        delta_x = find_x(U, find_z(L, minus_f));
        
        for (unsigned long long i = 0; i < n; ++i) 
            x_cur[i]= x_prev[i] + delta_x[i];
        iteration_cnt++;
        
        // calculate epsilon_k
        epsilon_k =  0;
        for (unsigned long long i = 0; i < n; ++i)
            epsilon_k = std::max(epsilon_k, std::abs(x_cur[i] - x_prev[i]));

        x_prev = x_cur;
    } while (epsilon_k > epsilon);

    return std::make_pair(x_cur[0], x_cur[1]);
}


int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open()) {

            double l1, r1, l2, r2;
            task_file >> l1 >> r1;
            task_file >> l2 >> r2;

            task_file >> epsilon;

            std::string test_number = argv[1];
            test_number.erase(0, 13); 
            std::string answer_file_name = "answer_" + test_number;
                        
            answer_file.open(answer_file_name);

            answer_file << "Fixed Point Iteration" << std::endl;
            auto [x1_fpi, x2_fpi] = fixed_point_iteration(l1, r1, l2, r2);
            answer_file << "x_1 = " << x1_fpi << std::endl;
            answer_file << "x_2 = " << x2_fpi << std::endl;
            answer_file << "Error: " << epsilon << std::endl;
            answer_file << "Iteration count: " << iteration_cnt << std::endl << std::endl;

            answer_file << "Newton's method (linearization)" << std::endl;
            auto [x1_nm, x2_nm] = newtons_method(r1, r2);
            answer_file << "x_1 = " << x1_nm << std::endl;
            answer_file << "x_2 = " << x2_nm << std::endl;
            answer_file << "Error: " << epsilon << std::endl;
            answer_file << "Iteration count: " << iteration_cnt;
            
            answer_file.close();

            std::cout << "Success! Answer was written to the file '" << answer_file_name << "'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
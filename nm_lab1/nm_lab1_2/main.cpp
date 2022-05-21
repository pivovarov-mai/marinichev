#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

unsigned long long n; // matrix size
std::ofstream answer_file_02;


void print_system(std::vector<std::vector<double>> matrix, std::vector<double> right_side) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_02 << matrix[i][j] << ' ';
        }
        answer_file_02 << right_side[i] << std::endl;
    }
}

void print_matrix(std::vector<std::vector<double>> matrix) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_02 << matrix[i][j] << ' ';
        }
        answer_file_02 << std::endl;
    }
}

void print_vector(std::vector<double> vector) {
    for (unsigned long long i = 0; i < n; ++i)
        answer_file_02 << 'x' << i+1 << " = " << vector[i] << std::endl;
}


// Tridiagonal matrix algorithm (Metod Progonki), aka Thomas algorithm (named after Llewellyn Thomas)
// |b_i| >= |a_i| + |c_i| -- reason why we can do this
std::vector<double> TDMA(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d) {
    std::vector<double> x(n);
    std::vector<double> P(n), Q(n);

    // moving up to find coefficients P_i and Q_i
    P[0] = - c[0] / b[0];
    Q[0] = d[0] / b[0];
    for (unsigned long long i = 1; i < n-1; ++i) {
        P[i] = - c[i] / (b[i] + a[i] * P[i-1]);
        Q[i] = (d[i] - a[i] * Q[i-1]) / (b[i] + a[i] * P[i-1]);
    } 
    x[n-1] = (d[n-1] - a[n-1]*Q[n-2]) / (b[n-1] + a[n-1]*P[n-2]);
    Q[n-1] = x[n-1];
    P[n-1] = 0;

    // moving down to find x_i
    for (unsigned long long i = n-1; i-- > 0;)
        x[i] = P[i]*x[i+1] + Q[i];

    return x;
}

int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream matrix_file(argv[1]);

        if (matrix_file.is_open()) {
            
            matrix_file >> n;

            std::vector<double> a(n); // digonal above
            std::vector<double> b(n); // main diagonal 
            std::vector<double> c(n); // digonal below
            std::vector<double> d(n); // right side
 
            // reading matrix 
            // first line
            a[0] = 0;
            matrix_file >> b[0] >> c[0];
            int zero; // variable to read zero elements of matrix
            for (unsigned long long i = 2; i < n; ++i) matrix_file >> zero;
            matrix_file >> d[0];

            // lines from second to (n-1)-th
            for (unsigned long long i = 1; i < n-1; ++i) {
                for (unsigned long long j = 0; j < i-1; ++j) matrix_file >> zero;
                for (unsigned long long j = i; j < i+1; ++j) matrix_file >> a[i] >> b[i] >> c[i];
                for (unsigned long long j = i+2; j < n; ++j) matrix_file >> zero;
                matrix_file >> d[i];
            }

            // n-th line
            c[n-1] = 0;
            for (unsigned long long i = 0; i < n-2; ++i) matrix_file >> zero;
            matrix_file >> a[n-1] >> b[n-1];
            matrix_file >> d[n-1];
            
            
            answer_file_02.open ("answer_02.txt");
            
            std::vector<double> X(n);
            X = TDMA(a, b, c, d);

            answer_file_02 << "X: " << std::endl;
            print_vector(X);

            answer_file_02.close();

            std::cout << "Success! Answer was written to the file 'answer_02.txt'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
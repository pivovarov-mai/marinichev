#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

unsigned long long n; // matrix size
unsigned long long permutation_cnt;
std::vector<unsigned long long> P(n); // vector to remember every permutation

std::ofstream answer_file_01;


void print_system(std::vector<std::vector<double>> matrix, std::vector<double> right_side) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_01 << matrix[i][j] << ' ';
        }
        answer_file_01 << right_side[i] << std::endl;
    }
}

void print_matrix(std::vector<std::vector<double>> matrix) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_01 << matrix[i][j] << ' ';
        }
        answer_file_01 << std::endl;
    }
}

void print_vector(std::vector<double> vector) {
    for (unsigned long long i = 0; i < n; ++i)
        answer_file_01 << 'x' << i+1 << " = " << vector[i] << std::endl;
}

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

void transpose(std::vector<std::vector<double>> &matrix) {
    std::vector<std::vector<double>> matrix_transposed(n);
    for (unsigned long long i = 0; i < n; ++i)
                matrix_transposed[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i)
        for (unsigned long long j = 0; j < n; ++j)
            matrix_transposed[i][j] = matrix[j][i];

    matrix = matrix_transposed;
}

void LU_decomposition(std::vector<std::vector<double>> a, std::vector<double> &B, std::vector<std::vector<double>> &L, std::vector<std::vector<double>> &U) {
    double mu;

    permutation_cnt = 0;
    P.resize(n);

    unsigned long long max_id;
    make_identity_matrix(L);
    for (unsigned long long k = 0; k < (n - 1); ++k) { // step
        max_id = find_max(a, k);
        P[k] = max_id; // remember permutation
        if (max_id != k) {
            swap_rows(a, k, max_id);
            std::swap(B[k], B[max_id]);
            swap_rows(L, k, max_id);
            swap_columns(L, k, max_id);
            permutation_cnt++;
        }

        for (unsigned long long i = (k + 1); i < n; ++i) { // index of a row
            mu = a[i][k] / a[k][k]; // coefficient to make k-th element of current row equal to zero
            L[i][k] = mu;
            for (unsigned long long j = k; j < n; ++j) { // index of a column
                a[i][j] -= mu * a[k][j]; // subtract k-th row multiplied by mu from all rows below it
            }
        }        
    }
    P[n-1] = n-1;
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

// x is A^-1
//  Ax = E
// LUx = E
//  Lz = E
//  Ux = z
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> l, std::vector<std::vector<double>> u, std::vector<std::vector<double>> &e) {
    make_identity_matrix(e);

    std::vector<std::vector<double>> z(n);
    std::vector<std::vector<double>> x(n);
    for (unsigned long long i = 0; i < n; ++i) {
        z[i].resize(n);
        x[i].resize(n);
    }

    for (unsigned long long i = 0; i < n; ++i) // iterate through eveyr column of identity matrix
        z[i] = find_z(l, e[i]);

    for (unsigned long long i = 0; i < n; ++i) // iterate through eveyr column of matrix z
        x[i] = find_x(u, z[i]);
    
    // undo all permutations with invertible matrix 
    transpose(x);
    for (unsigned long long k = n; k-- > 0; ) {
        if (k != P[k]) {            
            swap_columns(x, k, P[k]);
        }
    }

    return x;
}

double det(std::vector<std::vector<double>> U) {
    double det = 1;
    // multiplying the diagonal elements to get determinant
    for (unsigned long long i = 0; i < n; ++i)
        det *= U[i][i];

    if (permutation_cnt % 2 != 0)
        det *= (-1);
    return det;
}


int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream matrix_file(argv[1]);

        if (matrix_file.is_open()) {
            
            matrix_file >> n;

            std::vector<std::vector<double>> A(n);
            std::vector<double> B(n);
            std::vector<std::vector<double>> L(n);
            std::vector<std::vector<double>> U(n);
            std::vector<std::vector<double>> E(n);
            for (unsigned long long i = 0; i < n; ++i) {
                A[i].resize(n);
                L[i].resize(n);
                U[i].resize(n);
                E[i].resize(n);
            }

            // reading matrix 
            for (unsigned long long i = 0; i < n; ++i) {
                for (unsigned long long j = 0; j < n; ++j) {
                    matrix_file >> A[i][j];
                }
                matrix_file >> B[i];
            }
            
            LU_decomposition(A, B, L, U);

            
            answer_file_01.open ("answer_01.txt");

            answer_file_01 << "Matrix L: " << std::endl;
            print_matrix(L);
            answer_file_01 << std::endl;

            answer_file_01 << "Matrix U: " << std::endl;
            print_matrix(U);
            answer_file_01 << std::endl;
            
            std::vector<double> X(n);
            X = find_x(U, find_z(L, B));

            answer_file_01 << "X: " << std::endl;
            print_vector(X);
            answer_file_01 << std::endl;

            answer_file_01 << "A^-1: " << std::endl;
            print_matrix(inverse(L, U, E));
            answer_file_01 << std::endl;

            answer_file_01 << "det(A) = " << det(U);

            answer_file_01.close();

            std::cout << "Success! Answer was written to the file 'answer_01.txt'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
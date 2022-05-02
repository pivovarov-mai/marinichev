#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


unsigned long long n; // matrix size
unsigned long long iteration_cnt = 0;

std::ofstream answer_file_04;


void print_matrix(std::vector<std::vector<double>> matrix) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_04 << matrix[i][j] << ' ';
        }
        answer_file_04 << std::endl;
    }
}

void print_vectors(std::vector<std::vector<double>> vectors) {
    for (unsigned long long i = 0; i < n; ++i) {
        answer_file_04 << "V_" << i+1 << ": " << std::endl;
        for (unsigned long long j = 0; j < n; ++j)
            answer_file_04 << "v_" << i+1 << j+1 << " = " << vectors[j][i] << std::endl;
        answer_file_04 << std::endl;
    } 
}

void print_values(std::vector<double> values) {
    for (unsigned long long i = 0; i < n; ++i)
        answer_file_04 << "λ_" << i+1 << " = " << values[i] << std::endl;
}

void make_identity_matrix(std::vector<std::vector<double>> &I) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            if (i == j)
                I[i][j] = 1; 
            else
                I[i][j] = 0; 
        }
    }
}

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> matrix) {
    std::vector<std::vector<double>> matrix_transposed(n);
    for (unsigned long long i = 0; i < n; ++i)
        matrix_transposed[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i)
        for (unsigned long long j = 0; j < n; ++j)
            matrix_transposed[i][j] = matrix[j][i];

    return matrix_transposed;
}

std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> first_matrix, std::vector<std::vector<double>> second_matrix) {
    std::vector<std::vector<double>> product_matrix(n);
    for (unsigned long long i = 0; i < n; ++i)
        product_matrix[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            product_matrix[i][j] = 0;
            for (unsigned long long k = 0; k < n; ++k)
                product_matrix[i][j] += first_matrix[i][k] * second_matrix[k][j];
        }
    }
    return product_matrix;
}

std::vector<double> find_eigenvalues(std::vector<std::vector<double>> a) {
    std::vector<double> eigen_values(n);
    for (unsigned long long i = 0; i < n; ++i)
        eigen_values[i] = a[i][i];  

    return eigen_values;
}

// ||A||_2 = sqrt(sum_i,j;i<j A[i][j]^2)
double find_euclidean_norm(std::vector<std::vector<double>> a) {
    double matrix_norm = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        // only for non-diagonal elements
        for (unsigned long long j = 0; j < i; ++j)
            matrix_norm += a[i][j] * a[i][j];
    }
    matrix_norm = sqrt(matrix_norm);
    return matrix_norm;
}

std::pair<unsigned long long, unsigned long long> find_max(std::vector<std::vector<double>> a) {
    double max = std::fabs(a[1][0]);
    std::pair<unsigned long long, unsigned long long> max_id = {1, 0};
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < i; ++j) {
            if (std::fabs(a[i][j]) > max) {
                max = std::fabs(a[i][j]);
                max_id = {i, j};
            }
        }
    }
    return max_id; 
}

std::vector<double> jacobi_eigenvalue_algorithm(std::vector<std::vector<double>> a, std::vector<std::vector<double>> &v, double epsilon) {
	std::vector<std::vector<double>> u(n);
    std::vector<std::vector<double>> u_transposed;
    for (unsigned long long i = 0; i < n; ++i)
        u[i].resize(n);
    
    std::pair<unsigned long long, unsigned long long> cur_max_element_id;
    unsigned long long l, m;
    double phi, epsilon_k;
    bool do_initialize = true;

	do {
        cur_max_element_id = find_max(a);
        l = cur_max_element_id.first;
        m = cur_max_element_id.second;
        phi = 0.5 * atan((2 * a[l][m]) / (a[l][l] - a[m][m]));
        
        make_identity_matrix(u);
        u[l][l] = cos(phi);
        u[m][m] = cos(phi);
        u[l][m] = -sin(phi);
        u[m][l] = sin(phi);

        // A^k+1 = (U^k)^T * A^k * U^k
        u_transposed = transpose(u);
        a = multiply(u_transposed, a);
        a = multiply(a, u);

        // accumulate eigenvectors matrix
        if (do_initialize) {
            v = u;
            do_initialize = false;
        } else {
            v = multiply(v, u);
        }
   
        epsilon_k = find_euclidean_norm(a);

        iteration_cnt++;
	} while (epsilon_k > epsilon);

    return find_eigenvalues(a);
}

int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream matrix_file(argv[1]);

        if (matrix_file.is_open()) {
            
            matrix_file >> n;

            std::vector<std::vector<double>> A(n);
            std::vector<std::vector<double>> V(n);
            for (unsigned long long i = 0; i < n; ++i) {
                A[i].resize(n);
                V[i].resize(n);
            }

            // reading matrix 
            for (unsigned long long i = 0; i < n; ++i) {
                for (unsigned long long j = 0; j < n; ++j)
                    matrix_file >> A[i][j];
            }

            double error; // or epsilon
            matrix_file >> error;
                        
            answer_file_04.open ("answer_04.txt");

            std::vector<double> L(n);
            L = jacobi_eigenvalue_algorithm(A, V, error);

            answer_file_04 << "Eigenvalues: " << std::endl;
            print_values(L);
            answer_file_04 << std::endl;
           
            answer_file_04 << "Eigenvectors: " << std::endl;
            print_vectors(V);
            
            answer_file_04 << "Error: " << error << std::endl;

            answer_file_04 << "Iteration count: " << iteration_cnt;

            answer_file_04.close();

            std::cout << "Success! Answer was written to the file 'answer_04.txt'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
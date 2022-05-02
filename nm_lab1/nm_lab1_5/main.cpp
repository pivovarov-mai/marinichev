#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <cmath>


double epsilon; // or error
unsigned long long n; // matrix size
unsigned long long iteration_cnt = 0;

std::vector<std::vector<double>> A;
std::vector<std::complex<double>> eigenvalues;

std::ofstream answer_file_05;


void print_matrix(std::vector<std::vector<double>> matrix) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_05 << matrix[i][j] << ' ';
        }
        answer_file_05 << std::endl;
    }
}

std::string print_complex(std::complex<double> z) {
    if (abs(z.imag()) < epsilon)
        return std::to_string(z.real());
    else if (z.imag() < epsilon)
        return std::to_string(z.real()) + " - " + std::to_string(fabs(z.imag())) + "i";
    else
        return std::to_string(z.real()) + " + " + std::to_string(z.imag()) + "i";
}

void print_eigenvalues() {
    for (unsigned long long i = 0; i < n; ++i)
        answer_file_05 << "λ_" << i+1 << " = " << print_complex(eigenvalues[i]) << std::endl;
}


std::vector<std::vector<double>> E() {
    std::vector<std::vector<double>> E(n);
    for (unsigned long long i = 0; i < n; ++i)
        E[i].resize(n);
    
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            if (i == j)
                E[i][j] = 1; 
            else
                E[i][j] = 0; 
        }
    }
    
    return E;
}

void make_identity_matrix(std::vector<std::vector<double>> &I) {
    for (unsigned long long i = 0; i < n; ++i)
        I[i][i] = 1; 
}

std::vector<std::vector<double>> subtract(const std::vector<std::vector<double>> &first_matrix, const std::vector<std::vector<double>> &second_matrix) {
    std::vector<std::vector<double>> difference_matrix(n);
    for (unsigned long long i = 0; i < n; ++i)
        difference_matrix[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            difference_matrix[i][j] = first_matrix[i][j] - second_matrix[i][j];
        }
    }
    
    return difference_matrix;
}

std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>> &first_matrix, const std::vector<std::vector<double>> &second_matrix) {
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

std::vector<std::vector<double>> multiply_scalar(double scalar, const std::vector<std::vector<double>> &matrix) {
    std::vector<std::vector<double>> product_matrix(n);
    for (unsigned long long i = 0; i < n; ++i)
        product_matrix[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {            
            product_matrix[i][j] = scalar * matrix[i][j];
        }
    }
    
    return product_matrix;
}


std::pair<std::complex<double>, std::complex<double>> solve_characteristic_polynomial(double a_11, double a_12, double a_21, double a_22) {
    double a = 1.0;
    double b = -(a_11 + a_22);
    double c = a_11 * a_22 - a_12 * a_21;
    double discriminant_squared = b * b - 4.0 * a * c;
    if (discriminant_squared > epsilon) {
        std::complex<double> not_complex(NAN, NAN);
        return std::make_pair(not_complex, not_complex);
    }
    std::complex<double> discriminant(0.0, std::sqrt(-discriminant_squared));
    std::complex<double> x_1 = (-b + discriminant) / (2.0 * a);
    std::complex<double> x_2 = (-b - discriminant) / (2.0 * a);
    return std::make_pair(x_1, x_2);
}

void find_eigenvalues() {
    std::pair<std::complex<double>, std::complex<double>> complex_conjugate;
    for (unsigned long long i = 0; i < n; ++i) {
        if (i < n - 1 && (abs(A[i + 1][i]) > epsilon)) {
            complex_conjugate = solve_characteristic_polynomial(A[ i ][i], A[ i ][i+1], 
                                                                A[i+1][i], A[i+1][i+1]);
            if (std::isnan(complex_conjugate.first.real())) {
                ++i;
                continue;
            }
            eigenvalues[i] = complex_conjugate.first;
            eigenvalues[++i] = complex_conjugate.second;
        } else {
            eigenvalues[i] = A[i][i];
        }
    }
}


bool sub_sub_diagonal_converge() {
    double epsilon_k_2;
    for (unsigned long long i = 0; i < n; ++i) {
        epsilon_k_2 = 0;
        for (unsigned long long j = i + 2; j < n; ++j)
            epsilon_k_2 += A[j][i] * A[j][i];

        epsilon_k_2 = std::sqrt(epsilon_k_2);

        if (epsilon_k_2 > epsilon)
            return false;
    }
    return true;
}

bool converge() {
    if (!sub_sub_diagonal_converge())
        return false;

    std::vector<std::complex<double>> prev_eigenvalues(n);
    prev_eigenvalues = eigenvalues;
    find_eigenvalues();
    double epsilon_k_3;
    for (unsigned long long i = 0; i < n; ++i) {
        epsilon_k_3 = std::norm(eigenvalues[i] - prev_eigenvalues[i]);
        if (epsilon_k_3 > epsilon)
            return false;
    }
    return true;
}

double sign(double x) {
    if (x < epsilon)
        return -1.0;
    else if (x > epsilon)
        return 1.0;
    else
        return 0.0;
}

// ||v||_2 = sqrt(sum_i v[i]^2)
double find_euclidean_norm(const std::vector<double> &v) {
    double norm = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        norm += v[i] * v[i];
    }
    norm = std::sqrt(norm);
    return norm;
}

// v^T * v -> scalar
double multiply_vector_transposed_by_vector(const std::vector<double> &v) {
    double scalar = 0;
    for (unsigned long long i = 0; i < n; ++i)
        scalar += v[i] * v[i];
    
    return scalar;
}

// v * v^T -> matrix
std::vector<std::vector<double>> multiply_vector_by_vector_transposed(const std::vector<double> &v) {
    std::vector<std::vector<double>> matrix(n);
    for (unsigned long long i = 0; i < n; ++i)
        matrix[i].resize(n);
    
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            matrix[i][j] = v[i] * v[j];
        }
    }
    return matrix;
}

std::vector<std::vector<double>> make_householder_matrix(std::vector<double> v_k, unsigned long long id) {
    std::vector<double> v(n);
    v = v_k;
    v[id] += sign(v_k[id]) * find_euclidean_norm(v_k);
    
    return subtract(E(), multiply_scalar((2.0 / multiply_vector_transposed_by_vector(v)), multiply_vector_by_vector_transposed(v)));
}

void QR_decomposition(std::vector<std::vector<double>> &Q, std::vector<std::vector<double>> &R) {
    std::vector<std::vector<double>> H(n);
    for (unsigned long long i = 0; i < n; ++i) {
        H[i].resize(n);
    }

    make_identity_matrix(Q);
    R = A;
    for (unsigned long long i = 0; i < n - 1; ++i) {
        std::vector<double> v_k(n);
        for (unsigned long long j = i; j < n; ++j) {
            v_k[j] = R[j][i];
        }
        H = make_householder_matrix(v_k, i);
        Q = multiply(Q, H);
        R = multiply(H, R);
    }
}

void QR_algorithm() {
	do {
        std::vector<std::vector<double>> Q(n);
        std::vector<std::vector<double>> R(n);
        for (unsigned long long i = 0; i < n; ++i) {
            Q[i].resize(n);
            R[i].resize(n);
        }

        QR_decomposition(Q, R);
        A = multiply(R, Q);

        iteration_cnt++;
	} while (!converge());
}

int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream matrix_file(argv[1]);

        if (matrix_file.is_open()) {

            matrix_file >> n;

            A.resize(n);
            for (unsigned long long i = 0; i < n; ++i)
                A[i].resize(n);

            // reading matrix 
            for (unsigned long long i = 0; i < n; ++i) {
                for (unsigned long long j = 0; j < n; ++j)
                    matrix_file >> A[i][j];
            }
     
            matrix_file >> epsilon;

            std::string test_number = argv[1];
            test_number.erase(0, 13); 
            std::string answer_file_name = "answer_" + test_number;
                        
            answer_file_05.open(answer_file_name);

            eigenvalues.resize(n);
            QR_algorithm();

            answer_file_05 << "Eigenvalues: " << std::endl;
            print_eigenvalues();
            answer_file_05 << std::endl;
            
            answer_file_05 << "Error: " << epsilon << std::endl;

            answer_file_05 << "Iteration count: " << iteration_cnt;

            answer_file_05.close();

            std::cout << "Success! Answer was written to the file '" << answer_file_name << "'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
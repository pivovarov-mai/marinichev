#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


unsigned long long n; // matrix size
unsigned long long iteration_cnt;

std::ofstream answer_file_03;


void print_system(std::vector<std::vector<double>> matrix, std::vector<double> right_side) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_03 << matrix[i][j] << ' ';
        }
        answer_file_03 << right_side[i] << std::endl;
    }
}

void print_matrix(std::vector<std::vector<double>> matrix) {
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            answer_file_03 << matrix[i][j] << ' ';
        }
        answer_file_03 << std::endl;
    }
}

void print_vector(std::vector<double> vector) {
    for (unsigned long long i = 0; i < n; ++i)
        answer_file_03 << 'x' << i+1 << " = " << vector[i] << std::endl;
}

// alpha = E - D^-1 * A
std::vector<std::vector<double>> find_alpha_jacobi(std::vector<std::vector<double>> a) {
    std::vector<std::vector<double>> alpha(n);
    std::vector<std::vector<double>> d(n);
    for (unsigned long long i = 0; i < n; ++i) {
        alpha[i].resize(n);
        d[i].resize(n);
    }

    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            if (i == j) {
                d[i][j] = 1 / a[i][i];
                alpha[i][j] = 0;
            } else {
                alpha[i][j] = a[i][j];
            }
        }
    }

    double cur_sum = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            for (unsigned long long k = 0; k < n; ++k)
                cur_sum += d[i][k] * alpha[k][j];

            if (i != j) alpha[i][j] = -cur_sum;

            cur_sum = 0;
        }
    }
    return alpha;
}

// ||A||_inf = max_i (sum_j |A[i][j]|)
double find_matrix_norm(std::vector<std::vector<double>> a) {
    double matrix_norm = 0;
    double cur_sum = 0;
    for (unsigned long long j = 0; j < n; ++j)
        matrix_norm += fabs(a[0][j]);

    for (unsigned long long i = 1; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            cur_sum += fabs(a[i][j]);
        }
        if (cur_sum > matrix_norm) {
            matrix_norm = cur_sum;
        }
        cur_sum = 0;
    }
    return matrix_norm;
}

// x_k = alpha * x_k-1 + beta, where
// alpha = E - D^-1 * A
// beta = D^-1 * b
std::vector<double> fixed_point_iteration_jacobi(std::vector<std::vector<double>> a, std::vector<double> b, double epsilon) {
	std::vector<double> x_cur(n);  // x_k
    std::vector<double> x_prev(n); // x_k-1

    std::vector<std::vector<double>> alpha(n);
    for (unsigned long long i = 0; i < n; ++i)
        alpha[i].resize(n);
    alpha = find_alpha_jacobi(a);

    double alpha_norm;  // ||alpha||_inf = max_i (sum_j |alpha[i][j]|)
    alpha_norm = find_matrix_norm(alpha);
    
    // zero iteration
    iteration_cnt = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        x_prev[i] = b[i] / a[i][i];
    }

    double vector_norm; // ||x_k[i] - x_k-1[i]||_inf = max_i |x_k[i] - x_k-1[i]|
    double epsilon_k;
	do {
        // calculate x_k
		for (unsigned long long i = 0; i < n; ++i) {
			x_cur[i] = b[i];
			for (unsigned long long j = 0; j < n; ++j) {
				if (i != j)
					x_cur[i] -= a[i][j] * x_prev[j];
			}
			x_cur[i] /= a[i][i];
		}

        // calculate epsilon_k
        vector_norm = fabs(x_prev[0] - x_cur[0]);
		for (unsigned long long i = 0; i < n; ++i) {
			if (fabs(x_prev[i] - x_cur[i]) > vector_norm)
				vector_norm = fabs(x_prev[i] - x_cur[i]);
			x_prev[i] = x_cur[i];
		}
        epsilon_k =  (alpha_norm / (1 - alpha_norm)) * vector_norm;

        iteration_cnt++;
	} while (epsilon_k > epsilon);

    return x_cur;
}

std::vector<std::vector<double>> find_alpha(std::vector<std::vector<double>> a) {
    std::vector<std::vector<double>> alpha(n);
    for (unsigned long long i = 0; i < n; ++i)
        alpha[i].resize(n);

    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            if (i == j) {
                alpha[i][j] = 0;
            } else {
                alpha[i][j] = -(a[i][j] / a[i][i]);
            }
        }
    }
    return alpha;
}

std::vector<double> seidel_method(std::vector<std::vector<double>> a, std::vector<double> b, double epsilon) {
	std::vector<double> x_cur(n);  // x_k
    std::vector<double> x_prev(n); // x_k-1

    std::vector<std::vector<double>> alpha(n);
    for (unsigned long long i = 0; i < n; ++i)
        alpha[i].resize(n);
    alpha = find_alpha(a);

    double alpha_norm;  // ||alpha||_inf = max_i (sum_j |alpha[i][j]|)
    alpha_norm = find_matrix_norm(alpha);

    double c_norm;
    for (unsigned long long i = 0; i < n; ++i) {
        for (unsigned long long j = 0; j < n; ++j) {
            if (i <= j) 
                alpha[i][j] = 0;
        }
    }
    c_norm = find_matrix_norm(alpha);
    
    // zero iteration
    iteration_cnt = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        x_cur[i] = b[i] / a[i][i];
    }

    double vector_norm; // ||x_k[i] - x_k-1[i]||_inf = max_i |x_k[i] - x_k-1[i]|
    double epsilon_k;
    double sum_1, sum_2;
	do {
        // calculate x_k
        for (unsigned long long i = 0; i < n; ++i)
            x_prev[i] = x_cur[i];
        for (unsigned long long i = 0; i < n; ++i) {
            sum_1 = 0;
            for (unsigned long long j = 0; j < i; ++j)
                sum_1 += a[i][j] * x_cur[j];

            sum_2 = 0;  
            for (unsigned long long j = i+1; j < n; ++j)
                sum_2 += a[i][j] * x_prev[j];
   
            x_cur[i] = (b[i] - sum_1 - sum_2) / a[i][i];
        }

        // calculate epsilon_k
        vector_norm = fabs(x_prev[0] - x_cur[0]);
		for (unsigned long long i = 0; i < n; ++i) {
			if (fabs(x_prev[i] - x_cur[i]) > vector_norm)
				vector_norm = fabs(x_prev[i] - x_cur[i]);
		}
        epsilon_k =  (c_norm / (1 - alpha_norm)) * vector_norm;

        iteration_cnt++;
	} while (epsilon_k > epsilon);

    return x_cur;
}


int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream matrix_file(argv[1]);

        if (matrix_file.is_open()) {
            
            matrix_file >> n;

            std::vector<std::vector<double>> A(n);
            std::vector<double> B(n);
            for (unsigned long long i = 0; i < n; ++i)
                A[i].resize(n);

            // reading matrix 
            for (unsigned long long i = 0; i < n; ++i) {
                for (unsigned long long j = 0; j < n; ++j) {
                    matrix_file >> A[i][j];
                }
                matrix_file >> B[i];
            }

            double error; // or epsilon
            matrix_file >> error;
                        
            answer_file_03.open ("answer_03.txt");
           
            answer_file_03 << "Jacobi Method (Fixed Point Iteration)" << std::endl;

            std::vector<double> X(n);
            X = fixed_point_iteration_jacobi(A, B, error);

            answer_file_03 << "X: " << std::endl;
            print_vector(X);
            answer_file_03 << std::endl;

            answer_file_03 << "Error: " << error << std::endl;

            answer_file_03 << "Iteration count: " << iteration_cnt << std::endl;
            answer_file_03 << std::endl;

            answer_file_03 << "Seidel Method" << std::endl;

            X = seidel_method(A, B, error);

            answer_file_03 << "X: " << std::endl;
            print_vector(X);
            answer_file_03 << std::endl;

            answer_file_03 << "Error: " << error << std::endl;

            answer_file_03 << "Iteration count: " << iteration_cnt;

            answer_file_03.close();

            std::cout << "Success! Answer was written to the file 'answer_03.txt'" << std::endl;

        } else {
            std::cout << "Error: Unable to open file" << std::endl;
        } 
    } else {
        std::cout << "Syntax: ./solution [path to task file]" << std::endl;
    }

    return 0;
}
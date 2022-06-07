#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>


double epsilon; // or error
double a, b;    // segment [a, b] bounds
unsigned long long iteration_cnt = 0;

std::ofstream answer_file;


double log(double a, double b) {
    return log(b) / log(a);
}


double f(double x) {
    return std::pow(3, x) - 5 * x*x + 1;
}

double f_d(double x) {
    return std::pow(3, x) * std::log(3) - 10 * x;
}

double f_dd(double x) {
    return std::pow(3, x) * std::log(3) * std::log(3) - 10;
}


double phi_1(double x) {
    return (x > 0) ? std::sqrt((std::pow(3, x) + 1) / 5):
                    -std::sqrt((std::pow(3, x) + 1) / 5);
}

double phi_2(double x) {
    return log(3, 5 * x*x - 1);
}

double phi(double x) {
    if (phi_1(x) < 1)
        return phi_1(x);
    else
        return phi_2(x);
}


double phi_d_1(double x) {
    return (std::pow(3, x) * std::log(3)) / (2 * 5 * phi_1(x));             
}

double phi_d_2(double x) {
    return (10 * x) / ((5 * x*x - 1) * log(3));
}

double phi_d(double x) {
    if (phi_1(x) < 1)
        return phi_d_1(x);
    else
        return phi_d_2(x);
}


double fixed_point_iteration(double l, double r) {
    iteration_cnt = 0;
    double x_cur;
    double x_prev = (l + r) / 2;
    
    double q = std::max(std::abs(phi_d(l)), std::abs(phi_d(r)));
    double epsilon_k;
    do {
        x_cur = phi(x_prev);
        iteration_cnt++;
        epsilon_k = (q / (1 - q)) * std::abs(x_cur - x_prev);
        x_prev = x_cur;
    } while (epsilon_k > epsilon);

    return x_cur;
}

double newton_method(double l, double r) {
    iteration_cnt = 0;
    double x_cur;
    double x_prev;

    if (f(l) * f_dd(l) > 0)
        x_prev = l;
    else
        x_prev = r; 

    double epsilon_k;
    do {
        x_cur = x_prev - f(x_prev) / f_d(x_prev);
        iteration_cnt++;
        epsilon_k = std::abs(x_cur - x_prev);
        x_prev = x_cur;        
    } while (epsilon_k > epsilon);

    return x_cur;
}

int main(int argc, char *argv[]) {
    if (argc == 2) { // input from file
        std::ifstream task_file(argv[1]);

        if (task_file.is_open()) {

            task_file >> a >> b;
     
            task_file >> epsilon;

            std::string test_number = argv[1];
            test_number.erase(0, 13); 
            std::string answer_file_name = "answer_" + test_number;
                        
            answer_file.open(answer_file_name);

            answer_file << "Fixed Point Iteration" << std::endl;
            double x;
            x = fixed_point_iteration(a, b);
            answer_file << "x = " << x << std::endl;
            answer_file << "Error: " << epsilon << std::endl;
            answer_file << "Iteration count: " << iteration_cnt << std::endl << std::endl;

            answer_file << "Newton's method (linearization)" << std::endl;
            x = newton_method(a, b);
            answer_file << "x = " << x << std::endl;
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
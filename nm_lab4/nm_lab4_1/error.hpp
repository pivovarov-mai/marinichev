#ifndef ERROR_HPP
#define ERROR_HPP

#include <cmath>

const double EPS = 1e-9;

bool leq(double a, double b) { return (a < b) or (std::abs(b - a) < EPS); }

#endif /* ERROR_HPP */
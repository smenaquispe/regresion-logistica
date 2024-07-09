#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

double logistic_function(double beta0, double beta1, double x) {
    return 1.0 / (1.0 + exp(-(beta0 + beta1 * x)));
}

double cost_function(double beta0, double beta1, const std::vector<double>& x_values, const std::vector<double>& y_values) {
    int m = y_values.size();  // número de muestras
    double total_cost = 0.0;

    for (size_t i = 0; i < m; ++i) {
        double x = x_values[i];
        double y = y_values[i];
        double p = logistic_function(beta0, beta1, x);
        if (y == 1) {
            total_cost -= log(p);
        } else {
            total_cost -= log(1 - p);
        }
    }

    std::cout << total_cost << std::endl;

    return total_cost / m;
}

std::pair<double, double> gradient_descent(const std::vector<double>& x_values, const std::vector<double>& y_values, double beta0, double beta1, double learning_rate, int iterations) {
    int m = y_values.size();
    double current_cost;

    for (int i = 0; i < iterations; ++i) {
        double sum_errors_beta0 = 0.0;
        double sum_errors_beta1 = 0.0;

        for (size_t j = 0; j < m; ++j) {
            double x = x_values[j];
            double y = y_values[j];
            double p = logistic_function(beta0, beta1, x);
            double error = p - y;
            std::cout << "error " << error << std::endl;
            sum_errors_beta0 += error;
            sum_errors_beta1 += error * x;
        }

        beta0 -= learning_rate * sum_errors_beta0 / m;
        beta1 -= learning_rate * sum_errors_beta1 / m;

        current_cost = cost_function(beta0, beta1, x_values, y_values);
        std::cout << "Iteration " << i << ", Cost: " << current_cost << std::endl;
    }

    return std::make_pair(beta0, beta1);
}

int main() {
    double beta0 = 2.0;
    double beta1 = 1.0;
    std::vector<double> x_values = {1, 4, 6};  // valores de x
    std::vector<double> y_values = {1, 0, 0};  // valores de y (etiquetas)
    double learning_rate = 0.5;
    int iterations = 1;

    auto result = gradient_descent(x_values, y_values, beta0, beta1, learning_rate, iterations);
    std::cout << result.first << ", " << result.second << std::endl;

    auto test = logistic_function(result.first, result.second, 5);
    std::cout << test << std::endl;

    return 0;
}

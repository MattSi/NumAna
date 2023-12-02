#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "ch5_comm.h"

#define A (5.0)
static double fn1(double x, double y) {
    return A * y - A * x + 1;
}
static double fn1_solver(double x) {
    return (A * x - 1) / A;
}
static double fn1_orig(double x) {
    return std::exp(A * x) + x;
}

double rk4d(double(*fn)(double, double), double xn, double yn, double h)
{
    double k1 = h * fn(xn, yn);
    double k2 = h * fn(xn + 0.5 * h, yn + 0.5 * k1);
    double k3 = h * fn(xn + 0.5 * h, yn + 0.5 * k2);
    double k4 = h * fn(xn + h, yn + k3);

    double  yn_1 = yn + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    return yn_1;
}

double ab3d(double(*fn)(double, double), const std::vector<double> x, const std::vector<double> y, double h)
{
    if (y.size() < 3) {
        // There are no sufficient previous value to support this calculation
        return INFINITY;
    }
    int n = y.size() - 1;
    double yn = y.back();
    double yn_1 = yn + (h / 12.0) * (23 * fn(x[n], y[n]) - 16 * fn(x[n - 1], y[n - 1]) + 5 * fn(x[n - 2], y[n - 2]));
    
    return yn_1;
}

double am3d(double(*fn)(double, double), double(*solver)(double), const std::vector<double> x, const std::vector<double> y, double h)
{
    if (y.size() < 2) {
        return INFINITY;
    }
    int n = y.size() - 1;
    double xn = x.back();
    double yn = y.back();
    //double fn_1 = solver(xn+h);
    double yn_1 = yn + (h / 12.0) * (5.0 * (-A * (xn + h) + 1) + 8 * fn(xn, yn) - fn(x[n - 1], y[n - 1]));
    yn_1 /= (1 - 5 * h * A / 12);
    return yn_1;
}

void rk4d_test() {
    std::vector<double> xn_s;
    std::vector<double> yn_s;
    std::vector<double> yn_orig_s;

    double h = 0.01;
    double r_boundary = 1.0;
    xn_s.push_back(0.0);
    yn_s.push_back(1.0);
    yn_orig_s.push_back(1.0);

    while (true) {
        double last_x = xn_s.back();
        double last_y = yn_s.back();

        double curr_y = rk4d(fn1, last_x, last_y, h);
        double curr_x = last_x + h;
        if (curr_x < r_boundary) {
            xn_s.push_back(curr_x);
            yn_s.push_back(curr_y);

            yn_orig_s.push_back(fn1_orig(curr_x));
        }
        else {
            break;
        }
    }

    std::ofstream outputFile("ch5.txt");
    outputFile << std::fixed << std::setprecision(15);

    for (const auto& item : xn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_orig_s) {
        outputFile << item << " ";
    }

    outputFile << std::endl;


    outputFile.close();
}



void am3d_test() {
    std::vector<double> xn_s;
    std::vector<double> yn_s;
    std::vector<double> yn_orig_s;

    double h = 0.01;
    double r_boundary = 1.0;
    xn_s.push_back(0.0);
    yn_s.push_back(1.0);
    yn_orig_s.push_back(1.0);



    while (true) {
        double last_x = xn_s.back();
        double last_y = yn_s.back();

        double curr_y = rk4d(fn1, last_x, last_y, h);
        double curr_x = last_x + h;
        if (curr_x < r_boundary) {
            xn_s.push_back(curr_x);
            yn_s.push_back(curr_y);

            yn_orig_s.push_back(fn1_orig(curr_x));

            if (yn_s.size() == 3) {
                break;
            }
        }
        else {
            break;
        }
    }

    while (true) {
        double last_x = xn_s.back();
        double curr_y = ab3d(fn1, xn_s, yn_s, h);
        double curr_x = last_x + h;

        if (curr_x < r_boundary) {
            xn_s.push_back(curr_x);
            yn_s.push_back(curr_y);

            yn_orig_s.push_back(fn1_orig(curr_x));
        }
        else {
            break;
        }
    }

    std::ofstream outputFile("ch5_2.txt");
    outputFile << std::fixed << std::setprecision(15);

    for (const auto& item : xn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_orig_s) {
        outputFile << item << " ";
    }

    outputFile << std::endl;


    outputFile.close();
}

void ch5_driver() {
    // TODO Finish am3d testing 
    std::vector<double> xn_s;
    std::vector<double> yn_s;
    std::vector<double> yn_orig_s;

    double h = 0.01;
    double r_boundary = 1.0;
    xn_s.push_back(0.0);
    yn_s.push_back(1.0);
    yn_orig_s.push_back(1.0);



    while (true) {
        double last_x = xn_s.back();
        double last_y = yn_s.back();

        double curr_y = rk4d(fn1, last_x, last_y, h);
        double curr_x = last_x + h;
        if (curr_x < r_boundary) {
            xn_s.push_back(curr_x);
            yn_s.push_back(curr_y);

            yn_orig_s.push_back(fn1_orig(curr_x));

            if (yn_s.size() == 2) {
                break;
            }
        }
        else {
            break;
        }
    }

    while (true) {
        double last_x = xn_s.back();
        double curr_y = am3d(fn1, fn1_solver, xn_s, yn_s, h);
        double curr_x = last_x + h;

        if (curr_x < r_boundary) {
            xn_s.push_back(curr_x);
            yn_s.push_back(curr_y);

            yn_orig_s.push_back(fn1_orig(curr_x));
        }
        else {
            break;
        }
    }

    std::ofstream outputFile("ch5_3.txt");
    outputFile << std::fixed << std::setprecision(15);

    for (const auto& item : xn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_s) {
        outputFile << item << " ";
    }
    outputFile << std::endl;

    for (const auto& item : yn_orig_s) {
        outputFile << item << " ";
    }

    outputFile << std::endl;


    outputFile.close();
}
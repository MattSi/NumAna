#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include "ch3_comm.h"


double lagrange_interpolate(double* xi, double* yi, int n, double x)
{
    double ans = 0.0;
    for (int i = 0; i < n; i++) {
        if (std::abs(x - xi[i]) < 1e-15) {
            return yi[i];
        }
    }

    for (int i = 0; i < n; i++) {
        double denominator = 1.0;
        double numerator = 1.0;
        for (int j = 0; j < n; j++) {
            if (i == j)continue;
            numerator *= (x - xi[j]);
            denominator *= (xi[i] - xi[j]);
        }
        ans += (numerator * yi[i] / denominator);
    }
    return ans;
}

void cheybchev_poly_roots(int n, double* roots) {
    if (roots == NULL) {
        return;
    }
    for (int k = 0; k < n; k++) {
        roots[k] = std::cos((2 * k + 1) * M_PI / (2 * n));
    }
}

std::vector<double> divided_difference_table(const std::vector<double>& x, const std::vector<double> y)
{
    int n = x.size();
    std::vector<double> f = y;

    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--) {
            f[j] = (f[j] - f[j - 1]) / (x[j] - x[j - i]);
        }
    }
    return f;
}

double newton_interpolation(const std::vector<double>& x, const std::vector<double>& f, double target_x)
{
    double result = f[0];
    double term = 1.0;

    for (int i = 1; i < x.size(); i++) {
        term *= (target_x - x[i - 1]);
        result += term * f[i];
    }
    return result;
}

static double f1(double x) {
    return 1 / (25 + std::pow(x, 2));
}

static double f2(double x) {
    return x / (1 + std::pow(x, 4));
}

static double f3(double x) {
    return std::atan(x);
}


static void lagrange_test() {
    int n = 6;
    double* xi = (double*)std::calloc(n, sizeof(double));
    double* yi = (double*)std::calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        xi[i] = -1.0 + 2.0 * i / (n - 1);
        yi[i] = f1(xi[i]);
    }

    int m = 101;
    double* xioutput = (double*)std::calloc(m, sizeof(double));
    double* yioutput = (double*)std::calloc(m, sizeof(double));
    std::ofstream MyFile("ch3.csv");
    for (int i = 0; i < m; i++) {
        xioutput[i] = -1.0 + 2.0 * i / (m - 1);
        yioutput[i] = lagrange_interpolate(xi, yi, n, xioutput[i]);
        //fprintf(stdout, "%f\n", xioutput[i]);
        //MyFile << yioutput[i] << "\t";
    }
    //MyFile << std::endl;


    fprintf(stdout, "--------------------------------------------\n");
    std::memset(xi, 0, sizeof(double));
    std::memset(yi, 0, sizeof(double));

    n = 20;
    for (int i = 0; i < n; i++) {
        xi[i] = -5.0 + 10.0 * i / (n - 1);
        yi[i] = f3(xi[i]);
    }

    std::memset(xioutput, 0, sizeof(double));
    std::memset(yioutput, 0, sizeof(double));
    for (int i = 0; i < m; i++) {
        xioutput[i] = -5.0 + 10.0 * i / (m - 1);
        yioutput[i] = lagrange_interpolate(xi, yi, n, xioutput[i]);
        MyFile << xioutput[i] << "\t";
    }
    MyFile << std::endl;

    for (int i = 0; i < m; i++) {
        MyFile << yioutput[i] << "\t";
    }
    MyFile << std::endl;


    free(xi);
    free(yi);
    free(xioutput);
    free(yioutput);
    MyFile.close();
}


void cheybchev_test() {
    std::ofstream MyFile("ch3.csv");
    int n = 21;
    int m = 101;
    double* xi = (double*)std::calloc(n, sizeof(double));
    double* yi = (double*)std::calloc(n, sizeof(double));
    double* xioutput = (double*)std::calloc(m, sizeof(double));
    double* yioutput = (double*)std::calloc(m, sizeof(double));
    double* roots = (double*)calloc(n, sizeof(double));


    cheybchev_poly_roots(n, roots);
    for (size_t i = 0; i < n; i++)
    {
        fprintf(stdout, "%1.10f\n", roots[i]);
    }


    for (int i = 0; i < n; i++) {
        xi[i] = 5 * roots[i];
        yi[i] = f3(xi[i]);
    }


    for (int i = 0; i < m; i++) {
        xioutput[i] = -5.0 + 10.0 * i / (m - 1);
        yioutput[i] = lagrange_interpolate(xi, yi, n, xioutput[i]);
        MyFile << xioutput[i] << "\t";
    }
    MyFile << std::endl;

    for (int i = 0; i < m; i++) {
        MyFile << yioutput[i] << "\t";
    }
    MyFile << std::endl;


    free(xi);
    free(yi);
    free(xioutput);
    free(yioutput);
    MyFile.close();
    free(roots);
}
void ch3_driver() {
    
    double xi[5] = { 1,2,4,6,7 };
    double yi[5] = { 4,1,0,1,1 };

    std::vector<double> x = { 1,2,4,6,7 };
    std::vector<double> y = { 4,1,0,1,1 };
    auto f = divided_difference_table(x, y);
    for (const auto& c : f) std::cout << c << " ";
    std::cout << std::endl;

    std::cout << newton_interpolation(x, f, 3.5) << std::endl;
    
}


void ch3_driver2() {

}
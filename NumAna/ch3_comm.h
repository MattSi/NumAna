#pragma once
#include <iostream>
#include <vector>

double lagrange_interpolate(double *xi, double *yi, int n, double x);


void cheybchev_poly_roots(int n, double* roots);

std::vector<double> divided_difference_table(const std::vector<double>& x, const std::vector<double> y);

double newton_interpolation(const std::vector<double>& x, const std::vector<double>& f, double target_x);

void ch3_driver();

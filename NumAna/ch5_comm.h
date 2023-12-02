#pragma once
#include <vector>
/*
	ODE common solvers
	1) Rugne-Kutta 4 degree (Explicit)
	2) Adam-Bashforth method 4-degree (Explicit)
	3) Adam-Moulton method 4-degree (Implicit)
*/

double rk4d(double (*fn)(double, double), double xn, double yn, double h);

double ab3d(double (*fn)(double, double), const std::vector<double> x, const std::vector<double> y, double h);

double am3d(double (*fn)(double, double), double (*solver)(double),
	const std::vector<double> x, const std::vector<double> y, double h);
void ch5_driver();
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <memory.h>
#include "ch1_comm.h"

double f1(double x) {
	return std::pow(std::cos(x), 2) - x + 6;
}

double f2(double x) {
	return std::pow(x, 3) - 2;
}

double f3(double x) {
	return std::pow(x, 3) - x - 2;
}

double fn_1(double x) {
	return std::pow(x, 3) - 2 * x - 2;
}
double fn_1p(double x) {
	return std::pow(x, 2) * 3 - 2;
}


double bisect_info(double (*f)(double), double a, double b, double tor, double *mids, double *delta, int &N  ) {
	if (f(a) * f(b) > 0) {
		std::cerr << (stderr, "f(a)f(b) < 0 not satisfied!");
		return INFINITY;
	}
	
	int i = 0;
	double mid = 0.5 * (a + b);
	double fc = f(mid);
	double fa = f(a);
	double fb = f(b);
	mids[i] = mid;
	delta[i++] = fc;

	while (N > 0 && std::abs(fc) > tor) {

		if (fc * fa < 0) {
			b = mid;
			fb = fc;
		}
		else {
			a = mid;
			fa = fc;
		}
		mid = 0.5 * (a + b);
		fc = f(mid);
		mids[i] = mid;
		delta[i++] = fc;
		N--;
	}
	N = i;

	return mid;
}

double bisect(double (*f)(double), double a, double b, double tor, int N) {
	if (f(a) * f(b) > 0) {
		std::cerr<<(stderr, "f(a)f(b) < 0 not satisfied!");
		return INFINITY;
	}
	double mid = 0.5 * (a + b);
	double fc = f(mid);
	double fa = f(a);
	double fb = f(b);
	while (N > 0 && std::abs(fc) > tor) {
		
		if (fc * fa < 0) {
			b = mid;
			fb = fc;
		}
		else {
			a = mid;
			fa = fc;
		}
		mid = 0.5 * (a + b);
		fc = f(mid);
		N--;
	}

	return mid;
}


double fpi(double(*f)(double), double x0, double tor, double* mids, int& N)
{
	double xi = x0;
	double pi = 0;
	for (int i = 0; i < N; i++) {
		pi = f(xi);
		mids[i] = pi;
		if (std::abs(pi - xi) < tor) {
			N = i;
			break;
		}
		xi = pi;
	}
	return pi;
}

double newton_raphson(double(*f)(double), double(*fp)(double), double x0, double tor, double* mids, int& N)
{
	double xi = x0;
	double xi_1 = 0;
	double tmp_root = 0;
	for (int i = 0; i < N; i++) {
		
		double xi_1 = xi - f(xi) / fp(xi);
		mids[i] = xi_1;
		tmp_root = f(xi_1);
		if (std::abs(tmp_root) < tor) {
			// Find a root
			mids[++i] = tmp_root;
			N = i;
			break;
		}
		xi = xi_1;
	}
	return tmp_root;
}



double secant(double(*f)(double), double x0, double x1, double tor, double* mids, int& N)
{
	if (f(x0) * f(x1) > 0) {
		fprintf(stderr, "f(a)f(b) < 0 not satisfied!");
		return INFINITY;
	}
	double an = x0;
	double bn = x1;
	double tmp_root = 0;
	for (int i = 0; i < N; i++) {
		
		double pn_1 = bn - f(bn) * (bn - an) / (f(bn) - f(an));
		tmp_root = f(pn_1);
		mids[i] = pn_1;
		if (std::abs(tmp_root) < tor) {
			N = i;
			return pn_1;
		}
		if (f(an) * tmp_root < 0) {
			bn = pn_1;
		}
		else {
			an = pn_1;
		}
	}
	
	return an;
}

void ch1_driver() {
	int N = 1000;
	int n = N;
	double	*mids = (double*)calloc(sizeof(double), N), 
			*delta = (double*)calloc(sizeof(double), N);
	double ans_f1 = bisect_info(f1, 0, 8, 0.0000001, mids, delta, n);

	fprintf(stdout, "root of cos(x)^2-x+6=0 [0,8]:% 2.6f\n", ans_f1);
	fprintf(stdout, "========================================================\n");
	fprintf(stdout, "% 5s     %15s     % 17s\n", "NUMBER", "ROOT", "DELTA");
	for (int i = 0; i < n; i++) {
		fprintf(stdout, "% 5d    % 3.14f     % 3.14f\n", i+1, mids[i], delta[i]);
	}


	n = N;
	memset(mids, 0, N * sizeof(double));
	memset(delta, 0, N * sizeof(double));
	double ans_f2 = bisect_info(f2, 1, 2, 1e-8, mids, delta, n);
	fprintf(stdout, "\n\n\nroot of x^3-2=0 [1,2]:% 2.6f\n", ans_f2);
	fprintf(stdout, "========================================================\n");
	fprintf(stdout, "% 5s     %15s     % 17s\n", "NUMBER", "ROOT", "DELTA");
	for (int i = 0; i < n; i++) {
		fprintf(stdout, "% 5d    % 3.14f     % 3.14f\n", i + 1, mids[i], delta[i]);
	}

	//n = N;
	//memset(mids, 0, N * sizeof(double));
	//memset(delta, 0, N * sizeof(double));
	//double ans_f3 = fpi(f3, 1.4, 1e-8, mids, n);
	//fprintf(stdout, "\n\n\nroot of x^3-2*x-2=0 [1.4]:% 2.8f\n", ans_f3);
	//fprintf(stdout, "% 5s     %15s\n", "NUMBER", "ROOT");
	//for (int i = 0; i < n; i++) {
	//	fprintf(stdout, "% 5d    % 3.14f\n", i + 1, mids[i]);
	//}

	n = N;
	memset(mids, 0, N * sizeof(double));
	memset(delta, 0, N * sizeof(double));
	double ans_fn_1 = newton_raphson(fn_1, fn_1p, 5, 1e-9, mids, n);
	fprintf(stdout, "\n\n\nroot of x^3-2*x-2=0 [5]:% 2.8f\n", ans_fn_1);
	fprintf(stdout, "% 5s     %15s\n", "NUMBER", "ROOT");
	for (int i = 0; i < n; i++) {
		fprintf(stdout, "% 5d    % 3.14f\n", i + 1, mids[i]);
	}

	n = N;
	memset(mids, 0, N * sizeof(double));
	memset(delta, 0, N * sizeof(double));
	double ans_f_secant_1 = secant(fn_1, 1, 2, 1e-9, mids, n);
	fprintf(stdout, "\n\n\nroot of x^3-2*x-2=0 [1,2]:% 2.8f\n", ans_f_secant_1);
	fprintf(stdout, "% 5s     %15s\n", "NUMBER", "ROOT");
	for (int i = 0; i < n; i++) {
		fprintf(stdout, "% 5d    % 3.14f\n", i + 1, mids[i]);
	}
	free(mids);
	free(delta);
}
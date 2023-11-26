#pragma once
#include <math.h>

double bisect(double (*f)(double), double a, double b, double tor, int N=1000);
double bisect_info(double (*f)(double), double a, double b, double tor, double *mids, double *delta,int &N  );

// Fixed point iteration
double fpi(double (*f)(double), double x0, double tor, double* mids, int& N);

double newton_raphson(double (*f)(double), double (*fp)(double), double x0, double tor, double *mids, int &N);

double secant(double (*f)(double), double x0, double x1, double tor, double *mids, int &N);
void ch1_driver();
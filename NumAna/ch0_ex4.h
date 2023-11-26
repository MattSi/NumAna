#pragma once

double c0_4_1_orig(double x);
double c0_4_1_opt(double x);

double c0_4_1_2_orig(double x);
double c0_4_1_2_opt(double x);

double c0_4_3_opt(double a, double b);
void driver();

constexpr double INPUT[14] = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 ,1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, };
#include <stdio.h>
#include <math.h>
#include "ch0_ex4.h"

double c0_4_1_orig(double x) {
	return (1 - 1 / cos(x)) / pow(tan(x), 2);
}

double c0_4_1_opt(double x) {
	return -cos(x) / (cos(x) + 1);
}


double c0_4_1_2_orig(double x) {
    return (1 - pow((1 - x), 3)) / x;
}

double c0_4_1_2_opt(double x) {
    return 1 + (1 - x) * (2 - x);
}

double c0_4_3_opt(double a, double b) {
    return -pow(b, 2) / (a - sqrt(pow(a, 2) + pow(b, 2)));
}
void driver() {
    double y1, y2;
    for (int i = 0; i < sizeof(INPUT) / sizeof(INPUT[0]); i++) {
        double x = INPUT[i];
        y1 = c0_4_1_orig(x);
        y2 = c0_4_1_opt(x);
        fprintf(stdout, "% 3.14f    % 3.14f     % 3.14f\n", x, y1, y2);
    }

    fprintf(stdout, "===============================================================\n");
    for (int i = 0; i < sizeof(INPUT) / sizeof(INPUT[0]); i++) {
        double x = INPUT[i];
        y1 = c0_4_1_2_orig(x);
        y2 = c0_4_1_2_opt(x);
        fprintf(stdout, "% 3.14f    % 3.14f     % 3.14f\n", x, y1, y2);
    }

    fprintf(stdout, "===============================================================\n");
    fprintf(stdout, "%1.3e\n", c0_4_3_opt(-12345678987654321, 123));
}
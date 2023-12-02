#include <iostream>
#include <iomanip>
#include "ch2_comm.h"



VectorXd gauss_eliminate(MatrixXd A, VectorXd b)
{
    int current_col = 0;
    VectorXd tmp(A.rows(), 1);
    VectorXd root(A.rows(), 1);
    tmp.fill(0.0);
    root.fill(0.0);

    for (int i = 0; i < A.cols(); i++) {
        current_col = i;
        // 1. Find the largest A_ii to be the pivot
        int current_row = current_col;
        int pivot = A(current_row, current_col);
        int row_to_swap = -1;
        for (int p = current_row; p < A.rows(); p++) {
            if (A(p, current_col) > A(current_row, current_col)) {
                row_to_swap = p;
            }
        }

        if (row_to_swap != -1) {
            tmp = A.row(row_to_swap);
            A.row(row_to_swap) = A.row(current_row);
            A.row(current_row) = tmp;

            double tmp_b = b(row_to_swap, 0);
            b(row_to_swap, 0) = b(current_row, 0);
            b(current_row, 0) = tmp_b;
        }
        
        for (int j = current_row + 1; j < A.rows(); j++) {
            double factor = (-1) * A(j, i) / A(current_row, i);
            A.row(j) += A.row(current_row) * factor;
            b(j, 0) += b(current_row, 0) * factor;
        }
    }

    for (int i = A.rows() - 1; i >= 0; i--) {
        double tmp_sum = 0.0;
        for (int j = i + 1; j < A.cols(); j++) {
            tmp_sum += A(i, j) * root(j, 0);
        }
        root(i, 0) = (b(i, 0) - tmp_sum) / A(i, i);
    }
   
    return root;
}

void lu(MatrixXd A, MatrixXd& L, MatrixXd& U)
{
    MatrixXd middle_u = A;
    MatrixXd middle_l = MatrixXd::Identity(A.rows(), A.cols());
    for (int j = 0; j < A.cols() - 1; j++) {
        MatrixXd middle_j(A.rows(), A.cols());
        middle_j = MatrixXd::Identity(A.rows(), A.cols());

        for (int i = j + 1; i < A.rows(); i++) {
            middle_j(i, j) = -middle_u(i, j) / middle_u(j, j);
        }
        middle_u = middle_j * middle_u;
        middle_l += middle_j;
        
    }

    U = middle_u;
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                middle_l(i, j) = 1;
            }
            else {
                middle_l(i, j) = middle_l(i,j) * (-1);
            }
        }
        
    }
    L = middle_l;

    return ;
}

VectorXd jacobi_iteration(MatrixXd A, VectorXd b, bool disp_guess)
{
    Eigen::IOFormat VectorRowFmt(15, 0, ", ", " ", "", "", "[", "]");
    MatrixXd D = MatrixXd::Zero(A.rows(), A.cols());
    MatrixXd L = MatrixXd::Zero(A.rows(), A.cols());
    MatrixXd U = MatrixXd::Zero(A.rows(), A.cols());
    VectorXd guess_solution = VectorXd::Ones(b.rows());

    // Initialize Diagonal, Lower Triangular, Upper Triangular
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            if (i < j) {
                L(i, j) = A(i, j);
            }
            else if (i == j) {
                D(i, j) = A(i, j);
            }
            else {
                U(i, j) = A(i, j);
            }
        }
    }

    MatrixXd T = (-1) * D.inverse() * (L + U);
    MatrixXd C = D.inverse() * b;

    // First guess is ones;
    

    for (int i=1;;i++) {
        VectorXd last_guess = guess_solution;
        if (disp_guess) {
            std::cout << "ITERATION: " << i << " " << last_guess.format(VectorRowFmt) << std::endl;
        }
        guess_solution = T * last_guess + C;
        VectorXd delta = b - A * guess_solution;
        if (delta.norm() < 1e-9) {
            // We get the right answer
            if (disp_guess) {
                std::cout << "ITERATION: " << i + 1 << " " << guess_solution.format(VectorRowFmt) << std::endl;
            }
            break;
        }
    }

    return guess_solution;
}

VectorXd gauss_seidel(MatrixXd A, VectorXd b, bool disp_guess)
{
    Eigen::IOFormat VectorRowFmt(15, 0, ", ", " ", "", "", "[", "]");
    MatrixXd D = MatrixXd::Zero(A.rows(), A.cols());
    MatrixXd L = MatrixXd::Zero(A.rows(), A.cols());
    MatrixXd U = MatrixXd::Zero(A.rows(), A.cols());
    VectorXd guess_solution = VectorXd::Ones(b.rows());

    // Initialize Diagonal, Lower Triangular, Upper Triangular
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            if (i < j) {
                L(i, j) = A(i, j);
            }
            else if (i == j) {
                D(i, j) = A(i, j);
            }
            else {
                U(i, j) = A(i, j);
            }
        }
    }

    MatrixXd T = (-1) * (D + L).inverse() * U;
    MatrixXd C = (D + L).inverse() * b;

    // First guess is ones;
    for (int i = 1;; i++) {
        VectorXd last_guess = guess_solution;
        if (disp_guess) {
            std::cout << "ITERATION: " << i << " " << last_guess.format(VectorRowFmt) << std::endl;
        }
        
        guess_solution = T * last_guess + C;
        VectorXd delta = b - A * guess_solution;
        if (delta.norm() < 1e-9) {
            // We get the right answer
            if (disp_guess) {
                std::cout << "ITERATION: " << i+1 << " " << guess_solution.format(VectorRowFmt) << std::endl;
            }
            break;

        }
    }

    return guess_solution;
}

VectorXd back_substitution(MatrixXd A, VectorXd b) {
    VectorXd root(A.rows(), 1);
    root.fill(0.0);

    for (int i = A.rows() - 1; i >= 0; i--) {
        double tmp_sum = 0.0;
        for (int j = i + 1; j < A.cols(); j++) {
            tmp_sum += A(i, j) * root(j, 0);
        }
        root(i, 0) = (b(i, 0) - tmp_sum) / A(i, i);
    }

    return root;
}
VectorXd forward_substitution(MatrixXd A, VectorXd b) {
    VectorXd root(A.rows(), 1);
    root.fill(0.0);

    for (int i = 0; i < A.rows(); i++) {
        double tmp_sum = 0.0;
        for (int j = 0; j <= i - 1;j++) {
            tmp_sum += A(i, j) * root(j, 0);
        }
        root(i, 0) = (b(i, 0) - tmp_sum) / A(i, i);
    }

    return root;
}

# define nine 0.111111111111
# define eight 0.125
# define seven 0.142857142857
# define six 0.166666666666
# define five 0.2
# define four 0.25
# define three 0.333333333333
# define two 0.5
# define one 1.0

void ch2_driver()
{
    Eigen::IOFormat CleanFmt(9, 0, ", ", "\n", "[", "]");
    Eigen::IOFormat VectorRowFmt(15, 0, ", ", " ", "", "", "[", "]");
 
    std::string sep = "\n----------------------------------------\n";
    VectorXd b(5, 1);
    MatrixXd A(5, 5);

    b.fill(1);
    A.fill(0.0);
    A << nine, eight, seven, six, five,
        eight, seven, six, five, four,
        seven, six, five, four, three,
        six, five, four, three, two,
        five, four, three, two, one;

    VectorXd solution = gauss_eliminate(A, b);
    VectorXd delta = A * solution - b;
    std::cout << "Solution of equation Ax=b: (Gauss Elimination)\nA =\n" << A.format(CleanFmt) << "\nb =\n" << b.format(VectorRowFmt) << std::endl;
    std::cout << "Solution vector =\n" << solution.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta.norm() <<sep<< std::endl;

    VectorXd b2(4, 1);
    MatrixXd A2(4, 4);
    A2 << 7.2, 2.3, -4.4, 0.5,
        1.3, 6.3, -3.5, 2.8,
        5.6, 0.9, 8.1, -1.3,
        1.5, 0.4, 3.7, 5.9;
    b2 << 15.1, 1.8, 16.6, 36.9;
    VectorXd solution2 = gauss_eliminate(A2, b2);
    VectorXd delta2 = A2 * solution2 - b2;
    std::cout << "Solution of equation Ax=b: (Gauss Elimination)\nA =\n" << A2.format(CleanFmt) << "\nb =\n" << b2.format(VectorRowFmt) << std::endl;
    std::cout << "Solution vector =\n" << solution2.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta2.norm() << sep << std::endl;


    std::cout << "Solution  of equation Ax=b:(LU Decomposition)\nA =\n" << A.format(CleanFmt) << "\nb =\n" << b.format(VectorRowFmt) << std::endl;
    MatrixXd l = MatrixXd::Identity(5, 5);
    MatrixXd u = MatrixXd::Identity(5, 5);
    lu(A, l, u);
    VectorXd y = forward_substitution(l, b);
    VectorXd x = back_substitution(u, y);
    delta = A * x - b;
    std::cout << "L = \n" << l.format(CleanFmt) << std::endl;
    std::cout << "U = \n" << u.format(CleanFmt) << std::endl;
    std::cout << "Solution vector =\n" << x.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta.norm() << sep << std::endl;


    std::cout << "Solution  of equation Ax=b:(LU Decomposition)\nA =\n" << A2.format(CleanFmt) << "\nb =\n" << b2.format(VectorRowFmt) << std::endl;
    MatrixXd l2 = MatrixXd::Identity(4, 4);
    MatrixXd u2 = MatrixXd::Identity(4, 4);
    lu(A2, l2, u2);
    VectorXd y2 = forward_substitution(l2, b2);
    VectorXd x2 = back_substitution(u2, y2);
    delta2 = A2 * x2 - b2;
    std::cout << "L = \n" << l2.format(CleanFmt) << std::endl;
    std::cout << "U = \n" << u2.format(CleanFmt) << std::endl;
    std::cout << "Solution vector =\n" << x2.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta2.norm() << sep << std::endl;



    /*VectorXd b2(4, 1);
    MatrixXd A2(4, 4);
    A2 << 7.2, 2.3, -4.4, 0.5,
        1.3, 6.3, -3.5, 2.8,
        5.6, 0.9, 8.1, -1.3,
        1.5, 0.4, 3.7, 5.9;
    b2 << 15.1, 1.8, 16.6, 36.9;*/
    std::cout << "Solution of equation Ax=b: (Jacobi Iteration)\nA =\n" << A2.format(CleanFmt) << "\nb =\n" << b2.format(VectorRowFmt) << std::endl;
    solution2 = jacobi_iteration(A2, b2, true);
    delta2 = A2 * solution2 - b2;
    std::cout << "Solution vector =\n" << solution2.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta2.norm() << sep << std::endl;



    std::cout << "Solution of equation Ax=b: (Gauss_Seidel Iteration)\nA =\n" << A2.format(CleanFmt) << "\nb =\n" << b2.format(VectorRowFmt) << std::endl;
    solution2 = gauss_seidel(A2, b2, true);
    delta2 = A2 * solution2 - b2;
    std::cout << "Solution vector =\n" << solution2.format(VectorRowFmt) << std::endl;
    std::cout << "Delta: " << delta2.norm() << sep << std::endl;
}

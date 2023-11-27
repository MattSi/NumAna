#pragma once
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/*
	Suppose that the number of columns is equal to the number of rows

*/ 
VectorXd gauss_eliminate(MatrixXd A, VectorXd b);


/*
	Suppose that the number of columns is equal to rows
*/
void lu(MatrixXd A, MatrixXd& L, MatrixXd& U);



void ch2_driver();
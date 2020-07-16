
#include "lscm.h"
#include "vector_area_matrix.h"

#include <Eigen/SVD>
#include "eigs.h"
#include "repdiag.h"
#include "cotmatrix.h"
#include "massmatrix.h"
void lscm(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd & U)
{
	// Solve optimization as a generalized Eigen value problem
	//      min_U U'QU  subject to  U'BU = 1

	int n = V.rows(); 

	// Compute Q & B

	Eigen::SparseMatrix<double> A, L, Q;
	cgpl::cotmatrix(V, F, L);
	cgpl::repdiag(L, 2, Q);
	vector_area_matrix(F, A);
	Q = Q - A;

	Eigen::SparseMatrix<double> M, B;
	cgpl::massmatrix(V, F, M);
	cgpl::repdiag(M, 2, B);

	// Solve

	Eigen::MatrixXd sU;
	Eigen::VectorXd sS;
	cgpl::eigs(Q, B, 4, sU, sS);

	// Somehow, first 2 eigen value ~e^{-13}
	//      3, 4, ... looks more reasonable
	U.resize(n, 2);
	U << sU.col(2).topRows(n), sU.col(2).bottomRows(n);
}

#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "edge_lengths.h"
void cgpl::smooth(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & G,
	double lambda,
	Eigen::MatrixXd & U)
{

	Eigen::SparseMatrix<double> L;
	cotmatrix(V, F, L);

	Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
	massmatrix(V, F, M);

	// M * G = (M - lambda * L) * U 
	Eigen::SparseMatrix<double> A = -lambda * L;
	for(int i = 0; i < M.rows(); i++) {
		A.coeffRef(i, i) += M.diagonal()[i];
	}

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(A);
	U = solver.solve(M * G);
}

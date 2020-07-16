#include "arap_single_iteration.h"
void polar_svd3x3(const Eigen::Matrix3d &C,Eigen::Matrix3d &R)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix3d V = svd.matrixV();
	Eigen::Matrix3d U = svd.matrixU();
	R = U* V.transpose() ;
}
void arap_single_iteration(
	const Eigen::SparseMatrix<double> & Q,
	const Eigen::VectorXi & b,
	const Eigen::SparseMatrix<double> & K,
	const Eigen::MatrixXd & bc,
	Eigen::MatrixXd & U)
{
	Eigen::MatrixXd C = K.transpose() * U;
	Eigen::MatrixXd R(C.rows(), C.cols());

	for(int k = 0; k < Q.cols(); k++) {
		Eigen::Matrix3d C_k = C.block(k * 3, 0, 3, 3);
		Eigen::Matrix3d R_k;
		C_k.normalize();
		polar_svd3x3(C_k, R_k);
		R.block(k * 3, 0, 3, 3) = R_k;
	}

	Eigen::MatrixXd B= K * R;
	B = -B;

	U.resize(Q.cols(), 3);
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>solver;

	for (int i = 0; i < b.rows(); i++)
	{
		B(b(i),0) = bc(i, 0);
		B(b(i),1) = bc(i, 1);
		B(b(i),2) = bc(i, 2);
	}

	solver.compute(Q);
	Eigen::VectorXd X=solver.solve(B.col(0));
	Eigen::VectorXd Y=solver.solve(B.col(1));
	Eigen::VectorXd Z=solver.solve(B.col(2));

	U << X, Y, Z;
}

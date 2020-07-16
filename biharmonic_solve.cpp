#include "biharmonic_solve.h"
#include<Eigen/SparseLU>
void biharmonic_solve(
	const Eigen::SparseMatrix<double> & QQ,
	const Eigen::VectorXi & b,
	const Eigen::MatrixXd & bc,
	Eigen::MatrixXd & D)
{
	D.resize(QQ.cols(), 3);
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>solver;

	Eigen::VectorXd dx = Eigen::VectorXd::Zero(QQ.cols());
	Eigen::VectorXd dy = Eigen::VectorXd::Zero(QQ.cols());
	Eigen::VectorXd dz = Eigen::VectorXd::Zero(QQ.cols());
	for (int i = 0; i < b.rows(); i++)
	{
		dx(b(i)) = bc(i, 0);
		dy(b(i)) = bc(i, 1);
		dz(b(i)) = bc(i, 2);
	}

	solver.compute(QQ);
	Eigen::VectorXd Dx=solver.solve(dx);
	Eigen::VectorXd Dy=solver.solve(dy);
	Eigen::VectorXd Dz=solver.solve(dz);

	D << Dx, Dy, Dz;
}


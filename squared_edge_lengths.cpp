#include"squared_edge_lengths.h"
void cgpl::squared_edge_lengths(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::MatrixXd& L)
{
	L.resize(F.rows(), 3);

	int fn = F.rows();
	for (int i = 0; i < fn; i++)
	{
		L(i, 0) = (V.row(F(i, 1)) - V.row(F(i, 2))).squaredNorm();
		L(i, 1) = (V.row(F(i, 2)) - V.row(F(i, 0))).squaredNorm();
		L(i, 2) = (V.row(F(i, 0)) - V.row(F(i, 1))).squaredNorm();
	}
	L=L.array().eval();
}
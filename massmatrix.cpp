#include "massmatrix.h"
#include "doublearea.h"
#include "edge_lengths.h"
void cgpl::massmatrix(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
	Eigen::MatrixXd l;
	cgpl::edge_lengths(V, F, l);

	Eigen::VectorXd dblA(F.rows());
	cgpl::doublearea(l, dblA);

	Eigen::VectorXd areas = Eigen::VectorXd::Zero(F.maxCoeff() + 1);

	for (int i = 0; i < F.rows(); i++) 
	{
		Eigen::Vector3i f = F.row(i);
		double area = dblA[i];

		for (int j = 0; j < f.size(); j++) 
		{
			int v = f[j];
			areas[v] += area;
		}
	}

	M.resize(F.maxCoeff() + 1);
	M.diagonal() = 1 / (double)6 * areas;
}
void cgpl::massmatrix(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::SparseMatrix<double> & M)
{
	Eigen::DiagonalMatrix<double, Eigen::Dynamic>  MM;
	massmatrix(V, F, MM);
	int m = V.rows();

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	for (int i = 0; i < m; i++)
	{
		tripletList.push_back(T(i, i, MM.diagonal()[i]));
	}
	M.resize(m, m);
	M.setFromTriplets(tripletList.begin(), tripletList.end());
}

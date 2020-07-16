#include "biharmonic_precompute.h"

#include "cotmatrix.h"
#include "massmatrix.h"
void biharmonic_precompute(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & b,
	Eigen::SparseMatrix<double> & QQ)
{
	Eigen::SparseMatrix<double> L;
	cgpl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M;
	cgpl::massmatrix(V, F, M);

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(M);
	Eigen::SparseMatrix<double> M_inv_L = solver.solve(L).sparseView();
	Eigen::SparseMatrix<double> Q = L.transpose() * M_inv_L; // L^T M^-1 L

	typedef Eigen::Triplet<double>T;
	std::vector<T>tripleList;

	std::vector<bool> Known(V.rows(), false);
	for (int i = 0; i < b.rows(); i++)
	{
		Known[b(i)] = true;
		tripleList.push_back(T(b(i), b(i), 1));
	}


	for (int i = 0; i < Q.outerSize(); i++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(Q, i); it; ++it)
		{
			int row = it.row();
			int col = it.col();
			if (Known[row] == false)
			{
				tripleList.push_back(T(row, col, it.value()));
			}
		}
	}

	QQ.resize(V.rows(), V.rows());
	QQ.setFromTriplets(tripleList.begin(), tripleList.end());
}


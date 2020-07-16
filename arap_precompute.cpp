#include "arap_precompute.h"
#include "cotmatrix.h"
typedef Eigen::Triplet<double> T;
void arap_precompute(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & b,
	Eigen::SparseMatrix<double> & Q,
	Eigen::SparseMatrix<double> & K)
{
	Eigen::SparseMatrix<double> L;
	cgpl::cotmatrix(V, F, L);

	std::vector<T>tripleList;

	std::vector<bool> Known(V.rows(), false);
	for (int i = 0; i < b.rows(); i++)
	{
		Known[b(i)] = true;
		tripleList.push_back(T(b(i), b(i), 1));
	}


	for (int i = 0; i < L.outerSize(); i++)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
		{
			int row = it.row();
			int col = it.col();
			if (Known[row] == false)
			{
				tripleList.push_back(T(row, col, it.value()));
			}
		}
	}

	Q.resize(V.rows(), V.rows());
	Q.setFromTriplets(tripleList.begin(), tripleList.end());

	int nV = V.rows();
	std::vector<T> tripletList;

	for(int f = 0; f < F.rows(); f++) {
		for(int a = 0; a < 3; a++) {
			int i = F(f, (a + 1) % 3); 
			int j = F(f, (a + 2) % 3);

			Eigen::Vector3d diff = V.row(i) - V.row(j);
			Eigen::Vector3d eij = L.coeff(i, j) * diff / 6.0;
			for(int ki = 0; ki < 3; ki++) {
				int k = F(f, (a + ki) % 3);
				for(int B = 0; B < 3; B++) {
					tripletList.push_back(T(i, 3 * k + B, eij[B])); 
					tripletList.push_back(T(j, 3 * k + B, -eij[B]));
				}
			}
		}
	}

	K.resize(nV, 3 * nV);
	K.setFromTriplets(tripletList.begin(), tripletList.end());
}

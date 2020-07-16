#include "tutte.h"
#include "map_vertices_to_circle.h"
#include "boundary_loop.h"
#include "cotmatrix.h"

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
	Eigen::VectorXi boundary;
	cgpl::boundary_loop(F, boundary);

	Eigen::MatrixXd UV;
	cgpl::map_vertices_to_circle(V, boundary, UV);


	std::vector<bool> isBound(V.rows(), false);
	for (int i = 0; i < boundary.rows(); i++)
	{
		isBound[boundary(i, 0)] = true;
	}

	Eigen::SparseMatrix<double> L;
	// igl::cotmatrix(V, F, L);
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	for (int f = 0; f < F.rows(); f++) {
		for (int i = 0; i < 3; i++) {

			int ii = F(f, i);

			if (isBound[ii])
			{
				continue;
			}

			int j = (i + 1) % 3;
			int k = (i + 2) % 3;
			double norm_ij = (V.row(F(f, i)) - V.row(F(f, j))).norm();
			double norm_ik = (V.row(F(f, i)) - V.row(F(f, k))).norm();

			// case i != j
			tripletList.push_back(T(F(f, i), F(f, j), 1.0/norm_ij));
			tripletList.push_back(T(F(f, i), F(f, k), 1.0/norm_ik));
			// case i == j
			tripletList.push_back(T(F(f, i), F(f, i), -1.0/norm_ij));
			tripletList.push_back(T(F(f, i), F(f, i), -1.0/norm_ik));
		}
	}

	for (int i = 0; i < boundary.rows(); i++)
	{
		tripletList.push_back(T(boundary(i), boundary(i), 1));
	}

	L.resize(V.rows(), V.rows());
	L.setFromTriplets(tripletList.begin(), tripletList.end());

	U.resize(V.rows(), 2);
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>solver;

	Eigen::VectorXd Bu = Eigen::VectorXd::Zero(V.rows());
	Eigen::VectorXd Bv = Eigen::VectorXd::Zero(V.rows());
	for (int i = 0; i < boundary.rows(); i++)
	{
		Bu(boundary[i]) = UV(i, 0);
		Bv(boundary[i]) = UV(i, 1);
	}

	solver.compute(L);
	Eigen::VectorXd u=solver.solve(Bu);
	Eigen::VectorXd v=solver.solve(Bv);

	U << u, v;
}


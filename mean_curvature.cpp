#include "mean_curvature.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "per_vertex_normals.h"
void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
	Eigen::MatrixXd N;
	Eigen::SparseMatrix<double> L, M;

	cgpl::cotmatrix(V, F, L);
	cgpl::massmatrix(V, F, M);
	cgpl::per_vertex_normals(V, F, N);

	Eigen::MatrixXd Hn = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(M).solve(L * V);
	H = -(Hn * N.transpose()).diagonal();
}

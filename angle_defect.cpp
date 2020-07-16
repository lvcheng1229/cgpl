#include "angle_defect.h"
#include "internal_angles.h"

#include"squared_edge_lengths.h"
#include"massmatrix.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
	D = Eigen::VectorXd::Ones(V.rows()) * 2.0 * 3.14159265358979;
	
	Eigen::MatrixXd A;
	internal_angles(V,F, A);

	for(int i = 0; i < F.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			int v_i = F(i, j);
			D[v_i] -= A(i, j);
		}
	}

	Eigen::SparseMatrix<double> M;
	cgpl::massmatrix(V,F,M);
	Eigen::VectorXd Area = M.diagonal();

	D=D.array()/Area.array();
}

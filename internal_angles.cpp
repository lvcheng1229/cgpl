#include "internal_angles.h"
#include"squared_edge_lengths.h"
void internal_angles(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
  Eigen::MatrixXd & W)
{
	Eigen::MatrixXd l_sqr;
	cgpl::squared_edge_lengths(V, F, l_sqr);

	W.resize(l_sqr.rows(), l_sqr.cols());

	for(int i = 0; i < W.rows(); i++) {
		for(int j = 0; j < 3; j++) {
			double a2 = l_sqr(i, (j + 1) % 3);
			double b2 = l_sqr(i, (j + 2) % 3);
			double c2 = l_sqr(i, (j + 0) % 3);
			double a = sqrt(a2);
			double b = sqrt(b2);

			W(i, j) = acos((a2 + b2 - c2) / (2.0 * a * b));
		}
	}
}

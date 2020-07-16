
#ifndef CGPL_EDGE_LENGTHS_H
#define CGPL_EDGE_LENGTHS_H
#include <Eigen/Dense>

namespace cgpl
{
	void edge_lengths(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::MatrixXd& L);
}
#endif


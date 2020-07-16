#ifndef CGPL_SQUARED_EDGE_LENGTHS_H
#define CGPL_SQUARED_EDGE_LENGTHS_H
#include <Eigen/Dense>

namespace cgpl
{
	void squared_edge_lengths(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		Eigen::MatrixXd& L);
}

#endif
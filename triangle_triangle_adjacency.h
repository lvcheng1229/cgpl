#ifndef CGPL_TRIANGLE_TRIANGLE_ADJACENCY_H
#define CGPL_TRIANGLE_TRIANGLE_ADJACENCY_H
#include <Eigen/Core>
#include <vector>
namespace cgpl
{
	void triangle_triangle_adjacency(
		const Eigen::MatrixXi& F,
		Eigen::MatrixXi& TT);
}
#endif
#ifndef CGPL_BOUNDARY_LOOP_H
#define CGPL_BOUNDARY_LOOP_H

#include <Eigen/Dense>
#include <vector>

namespace cgpl
{
	// Compute list of ordered boundary loops for a manifold mesh.
	//
	// Inputs:
	//   F  #V by dim list of mesh faces
	// Outputs:
	//   L  list of loops where L[i] = ordered list of boundary vertices in loop i
	//
	void boundary_loop(
		const Eigen::MatrixXi & F,
		std::vector<std::vector<int> >& L);

	//找边界顶点数最大的一个Loop
	void boundary_loop(
		const Eigen::MatrixXi & F,
		std::vector<int> & L);

	void boundary_loop(
		const Eigen::MatrixXi & F,
		Eigen::VectorXi & L);
}
#endif
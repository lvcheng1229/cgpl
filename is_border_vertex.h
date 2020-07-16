#ifndef CGPL_IS_BORDER_VERTEX_H
#define CGPL_IS_BORDER_VERTEX_H

#include <Eigen/Core>
#include <vector>
namespace cgpl
{
	// Determine vertices on open boundary of a (manifold) mesh with triangle
	// faces F
	//
	// Inputs:
	//   V  #V by dim list of vertex positions 
	//   F  #F by 3 list of triangle indices
	// Returns #V vector of bools revealing whether vertices are on boundary
	//
	// Known Bugs: - assumes mesh is edge manifold
	// 
	std::vector<bool> is_border_vertex(
		const Eigen::MatrixXi &F);
}
#endif
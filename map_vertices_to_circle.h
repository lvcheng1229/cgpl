#ifndef CGPL_MAP_VERTICES_TO_CIRCLE_H
#define CGPL_MAP_VERTICES_TO_CIRCLE_H
#include "PI.h"

#include <Eigen/Dense>
#include <vector>

namespace cgpl
{

	// Map the vertices whose indices are in a given boundary loop (bnd) on the
	// unit circle with spacing proportional to the original boundary edge
	// lengths.
	//
	// Inputs:
	//   V  #V by dim list of mesh vertex positions
	//   b  #W list of vertex ids
	// Outputs:
	//   UV   #W by 2 list of 2D position on the unit circle for the vertices in b
	void map_vertices_to_circle(
		const Eigen::MatrixXd& V,
		const Eigen::VectorXi& bnd,
		Eigen::MatrixXd& UV);
}
#endif

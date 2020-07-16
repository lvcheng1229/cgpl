#ifndef CGPL_PER_VERTEX_NORMALS_H
#define CGPL_PER_VERTEX_NORMALS_H

#include <Eigen/Core>

namespace cgpl
{
	enum PerVertexNormalsWeightingType
	{
		// Incident face normals have uniform influence on vertex normal
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM = 0,
		// Incident face normals are averaged weighted by area
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA = 1,
		// Incident face normals are averaged weighted by incident angle of vertex
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE = 2,
		// Area weights
		PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT = 3,
		NUM_PER_VERTEX_NORMALS_WEIGHTING_TYPE = 4
	};

	void per_vertex_normals(  
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const PerVertexNormalsWeightingType weighting,
		Eigen::MatrixXd &N);

	void per_vertex_normals(  
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::MatrixXd &N);
}

#endif

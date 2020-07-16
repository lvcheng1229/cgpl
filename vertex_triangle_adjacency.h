#include <Eigen/Dense>
#include <vector>
namespace cgpl
{
	//[n] number of vertex
	//[F] matrix of triangle face
	//[VF] adjacency triangle list of vertex V
	void vertex_triangle_adjacency(int n,
		const Eigen::MatrixXi &F,
		std::vector<std::vector<int> > &VF,
		std::vector<std::vector<int> > &VFi);

	//[Ni] number of adjacency triangle
	void vertex_triangle_adjacency(int n,
		const Eigen::MatrixXi &F,
		std::vector<std::vector<int> > &VF,
		std::vector<int> &Ni);
}
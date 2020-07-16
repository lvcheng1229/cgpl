#include "is_border_vertex.h"
#include"triangle_triangle_adjacency.h"
std::vector<bool> cgpl::is_border_vertex(
	const Eigen::MatrixXi &F)
{
	Eigen::MatrixXi FF;
	cgpl::triangle_triangle_adjacency(F,FF);
	std::vector<bool> ret(F.maxCoeff()+1);
	for(unsigned i=0; i<ret.size();++i)
		ret[i] = false;

	for(unsigned i=0; i<F.rows();++i)
		for(unsigned j=0;j<F.cols();++j)
			if(FF(i,j) == -1)
			{
				ret[F(i,j)]       = true;
				ret[F(i,(j+1)%F.cols())] = true;
			}
	return ret;
}
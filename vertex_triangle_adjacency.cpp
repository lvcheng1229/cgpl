#include"vertex_triangle_adjacency.h"
void cgpl::vertex_triangle_adjacency(int n,
	const Eigen::MatrixXi &F,
	std::vector<std::vector<int> > &VF,
	std::vector<std::vector<int> > &VFi)
{
	VF.clear();
	VFi.clear();

	VF.resize(n);
	VFi.resize(n);


	for(int fi=0; fi<F.rows(); ++fi)
	{
		for(int i = 0; i < F.cols(); ++i)
		{
			VF[F(fi,i)].push_back(fi);
			VFi[F(fi,i)].push_back(i);
		}
	}
}

void cgpl::vertex_triangle_adjacency(int n,
	const Eigen::MatrixXi &F,
	std::vector<std::vector<int> > &VF,
	std::vector<int> &Ni)
{
	VF.clear();
	Ni.clear();

	VF.resize(n);
	Ni.resize(n);

	for (auto t : Ni)
	{
		t = 0;
	}

	for(int fi=0; fi<F.rows(); ++fi)
	{
		for(int i = 0; i < F.cols(); ++i)
		{
			VF[F(fi,i)].push_back(fi);
			Ni[F(fi,i)]++;
		}
	}
}

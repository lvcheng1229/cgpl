#include"triangle_triangle_adjacency.h"
#include"vertex_triangle_adjacency.h"
void cgpl::triangle_triangle_adjacency(
	const Eigen::MatrixXi& F,
	Eigen::MatrixXi& TT)
{
	const int n = F.maxCoeff()+1;

	std::vector<int>Ni;
	std::vector<std::vector<int> > VF;
	cgpl::vertex_triangle_adjacency(n, F, VF, Ni);

	int fn = F.rows();
	TT = Eigen::MatrixXi::Constant(fn, 3, -1);
	// Loop over faces
	for (int i = 0; i < fn; i++)
	{
		// Loop over corners
		for (int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3);
			// Loop over face neighbors incident on this corner
			for (int k = 0; k < Ni[a]; k++)
			{
				int ff = VF[a][k];
				// Not this face
				if (ff != i)
				{
					// Face neighbor also has [vi,vin] edge
					if (F(ff,0) == b || F(ff,1) == b || F(ff,2) == b)
					{
						TT(i,j) = ff;
						break;
					}
				}
			}
		}
	}
	
}
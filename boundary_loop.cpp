#include "boundary_loop.h"
#include "triangle_triangle_adjacency.h"
#include "vertex_triangle_adjacency.h"
#include "is_border_vertex.h"
#include <set>
void cgpl::boundary_loop(
	const Eigen::MatrixXi & F,
	std::vector<std::vector<int> >& L)
{
	using namespace std;
	using namespace Eigen;

	if(F.rows() == 0)
		return;

	MatrixXi TT;
	vector<std::vector<int>> VF, VFi;
	triangle_triangle_adjacency(F,TT);
	vertex_triangle_adjacency(F.maxCoeff()+1,F,VF,VFi);

	vector<bool> unvisited = is_border_vertex(F);
	set<int> unseen;

	for (size_t i = 0; i < unvisited.size(); ++i)
	{
		if (unvisited[i])
			unseen.insert(unseen.end(),i);
	}

	while (!unseen.empty())
	{
		vector<int> l;

		// Get first vertex of loop
		int start = *unseen.begin();
		unseen.erase(unseen.begin());
		unvisited[start] = false;
		l.push_back(start);

		bool done = false;
		while (!done)
		{
			// Find next vertex
			bool newBndEdge = false;
			int v = l[l.size()-1];
			int next;
			for (int i = 0; i < (int)VF[v].size() && !newBndEdge; i++)
			{
				int fid = VF[v][i];

				if (TT.row(fid).minCoeff() < 0.) // Face contains boundary edge
				{
					int vLoc = -1;
					if (F(fid,0) == v) vLoc = 0;
					if (F(fid,1) == v) vLoc = 1;
					if (F(fid,2) == v) vLoc = 2;

					int vNext = F(fid,(vLoc + 1) % F.cols());

					newBndEdge = false;
					if (unvisited[vNext] && TT(fid,vLoc) < 0)
					{
						next = vNext;
						newBndEdge = true;
					}
				}
			}

			if (newBndEdge)
			{
				l.push_back(next);
				unseen.erase(next);
				unvisited[next] = false;
			}
			else
				done = true;
		}
		L.push_back(l);
	}

}

void cgpl::boundary_loop(
	const Eigen::MatrixXi & F,
	std::vector<int> & L)
{
	using namespace Eigen;
	using namespace std;

	if(F.rows() == 0)
		return;

	vector<vector<int> > Lall;
	boundary_loop(F,Lall);

	int idxMax = -1;
	size_t maxLen = 0;
	for (size_t i = 0; i < Lall.size(); ++i)
	{
		if (Lall[i].size() > maxLen)
		{
			maxLen = Lall[i].size();
			idxMax = i;
		}
	}

	//Check for meshes without boundary
	if (idxMax == -1)
	{
		L.clear();
		return;
	}

	L.resize(Lall[idxMax].size());
	for (size_t i = 0; i < Lall[idxMax].size(); ++i)
	{
		L[i] = Lall[idxMax][i];
	}
}

void cgpl::boundary_loop(
	const Eigen::MatrixXi & F,
	Eigen::VectorXi & L)
{
	using namespace Eigen;
	using namespace std;

	if(F.rows() == 0)
		return;

	vector<int> Lvec;
	boundary_loop(F,Lvec);

	L.resize(Lvec.size());
	for (size_t i = 0; i < Lvec.size(); ++i)
		L(i) = Lvec[i];
}
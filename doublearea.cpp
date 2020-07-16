#include "doublearea.h"
#include"edge_lengths.h"
#include<math.h>
void cgpl::doublearea(
	const Eigen::MatrixXd & l,
	Eigen::VectorXd & dblA)
{
	int m = l.rows();

	dblA.resize(l.rows(),1);

	for (int i = 0; i < m; i++)
	{
		const double arg =
			(l(i,0)+(l(i,1)+l(i,2)))*
			(l(i,2)-(l(i,0)-l(i,1)))*
			(l(i,2)+(l(i,0)-l(i,1)))*
			(l(i,0)+(l(i,1)-l(i,2)));
		dblA(i) = 2.0*0.25*sqrt(arg);
	}
}

void cgpl::doublearea(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::VectorXd & A)
{
	Eigen::MatrixXd L;
	edge_lengths(V, F, L);

	doublearea(L, A);
}
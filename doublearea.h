#ifndef CGPL_DOUBLEAREA_H
#define CGPL_DOUBLEAREA_H
#include <Eigen/Dense>
namespace cgpl
{
	void doublearea(
		const Eigen::MatrixXd & l,
		Eigen::VectorXd & dblA);

	void doublearea(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::VectorXd & A);
}
#endif
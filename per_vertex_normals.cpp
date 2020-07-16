#include"per_vertex_normals.h"
#include"internal_angles.h"
#include"doublearea.h"
void cgpl::per_vertex_normals(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const PerVertexNormalsWeightingType weighting,
	Eigen::MatrixXd &N)
{
	Eigen::MatrixXd FN;

	FN.resize(F.rows(),3);
	// loop over faces
	int Frows = F.rows();

	Eigen::Matrix<double, 1, 3>Z(0, 0, 0);

	for(int i = 0; i < Frows;i++)
	{
		const Eigen::Matrix<double, 1, 3> v1 = V.row(F(i,1)) - V.row(F(i,0));
		const Eigen::Matrix<double, 1, 3> v2 = V.row(F(i,2)) - V.row(F(i,0));

		FN.row(i) = v1.cross(v2);//.normalized();
		double r = FN.row(i).norm();
		if(r == 0)
		{
			FN.row(i) = Z;
		}else
		{
			FN.row(i) /= r;
		}
	}

	using namespace std;
	// Resize for output
	N.setZero(V.rows(),3);

	Eigen::MatrixXd W;

	W.resize(F.rows(), 3);
	switch(weighting)
	{
	case PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM:
		W.setConstant(1.);
		break;
	default:
		assert(false && "Unknown weighting type");
	case PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT:
	case PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA:
	{
		Eigen::VectorXd A;
		doublearea(V,F,A);
		W = A.replicate(1,3);
		break;
	}
	case PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE:
		internal_angles(V,F,W);
		break;
	}

	// loop over faces
	for(int i = 0;i<F.rows();i++)
	{
		// throw normal at each corner
		for(int j = 0; j < 3;j++)
		{
			N.row(F(i,j)) += W(i,j) * FN.row(i);
		}
	}

	// take average via normalization
	N.rowwise().normalize();
}

void cgpl::per_vertex_normals(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	Eigen::MatrixXd &N)
{
	per_vertex_normals(V, F, PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT, N);
}

#ifndef CGPL_EIGS_H
#define CGPL_EIGS_H

#include <Eigen/Core>
#include <Eigen/Sparse>


#include "sort.h"
#include <iostream>
namespace cgpl
{
	// Act like MATLAB's eigs function. Compute the first/last k eigen pairs of
	// the generalized eigen value problem:
	//
	//     A u = s B u
	//
	// Solutions are approximate and sorted. 
	//
	// Ideally one should use ARPACK and the Eigen unsupported ARPACK module.
	// This implementation does simple, naive power iterations.
	//
	// Inputs:
	//   A  #A by #A symmetric matrix
	//   B  #A by #A symmetric positive-definite matrix
	//   k  number of eigen pairs to compute
	//   type  whether to extract from the high or low end
	// Outputs:
	//   sU  #A by k list of sorted eigen vectors (descending)
	//   sS  k list of sorted eigen values (descending)
	//
	// Known issues:
	//   - only the 'sm' small magnitude eigen values are well supported
	//   
	template <
		typename Atype,
		typename Btype,
		typename DerivedU,
		typename DerivedS>
		bool eigs(
			const Eigen::SparseMatrix<Atype> & A,
			const Eigen::SparseMatrix<Btype> & iB,
			const size_t k,
			Eigen::PlainObjectBase<DerivedU> & sU,
			Eigen::PlainObjectBase<DerivedS> & sS);
}

template <
	typename Atype,
	typename Btype,
	typename DerivedU,
	typename DerivedS>
	bool cgpl::eigs(
		const Eigen::SparseMatrix<Atype> & A,
		const Eigen::SparseMatrix<Btype> & iB,
		const size_t k,
		Eigen::PlainObjectBase<DerivedU> & sU,
		Eigen::PlainObjectBase<DerivedS> & sS)
{
	using namespace Eigen;
	using namespace std;
	const size_t n = A.rows();
	assert(A.cols() == n && "A should be square.");
	assert(iB.rows() == n && "B should be match A's dims.");
	assert(iB.cols() == n && "B should be square.");
	DerivedU U(n,k);
	DerivedS S(k,1);
	typedef Atype Scalar;
	typedef Eigen::Matrix<typename DerivedU::Scalar,DerivedU::RowsAtCompileTime,1> VectorXS;
	// Rescale B for better numerics
	const Scalar rescale = std::abs(iB.diagonal().maxCoeff());
	const Eigen::SparseMatrix<Btype> B = iB/rescale;

	Scalar tol = 1e-4;
	Scalar conv = 1e-14;
	int max_iter = 100;
	int i = 0;
	while(true)
	{
		// Random initial guess
		VectorXS y = VectorXS::Random(n,1);
		Scalar eff_sigma = 0;
		if(i>0)
		{
			eff_sigma = 1e-8+std::abs(S(i-1));
		}
		// whether to use rayleigh quotient method
		bool ray = false;
		Scalar err = std::numeric_limits<Scalar>::infinity();
		int iter;
		Scalar sigma = std::numeric_limits<Scalar>::infinity();
		VectorXS x;
		for(iter = 0;iter<max_iter;iter++)
		{
			if(i>0 && !ray)
			{
				// project-out existing modes
				for(int j = 0;j<i;j++)
				{
					const VectorXS u = U.col(j);
					y = (y - u*u.dot(B*y)/u.dot(B * u)).eval();
				}
			}
			// normalize
			x = y/sqrt(y.dot(B*y));

			// current guess at eigen value
			sigma = x.dot(A*x)/x.dot(B*x);
			//x *= sigma>0?1.:-1.;

			Scalar err_prev = err;
			err = (A*x-sigma*B*x).array().abs().maxCoeff();
			if(err<conv)
			{
				break;
			}
			if(ray || err<tol)
			{
				eff_sigma = sigma;
				ray = true;
			}

			Scalar tikhonov = std::abs(eff_sigma)<1e-12?1e-10:0;
			SimplicialLDLT<SparseMatrix<Scalar> > solver;
			const SparseMatrix<Scalar> C = A-eff_sigma*B+tikhonov*B;
			solver.compute(C);
			switch(solver.info())
			{
			case Eigen::Success:
				break;
			case Eigen::NumericalIssue:
				cerr<<"Error: Numerical issue."<<endl;
				return false;
			default:
				cerr<<"Error: Other."<<endl;
				return false;
			}
			const VectorXS rhs = B*x;
			y = solver.solve(rhs);
		}
		if(iter == max_iter)
		{
			cerr<<"Failed to converge."<<endl;
			return false;
		}
		if(
			i==0 || 
			(S.head(i).array()-sigma).abs().maxCoeff()>1e-14 ||
			((U.leftCols(i).transpose()*B*x).array().abs()<=1e-7).all()
			)
		{
			U.col(i) = x;
			S(i) = sigma;
			i++;
			if(i == k)
			{
				break;
			}
		}
	}
	// finally sort
	VectorXi I;
	cgpl::sort(S,false,sS,I);

	sU.resize(U.rows(), U.cols());
	for (int i = 0; i < U.cols(); i++)
	{
		sU.col(i) = U.col(I(i));
	}
	sS /= rescale;
	sU /= sqrt(rescale);
	return true;
}

#endif

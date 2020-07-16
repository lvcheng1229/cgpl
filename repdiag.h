// This file is part of libigl, a simple c++ geometry processing library.
// 
#ifndef CGPL_REPDIAG_H
#define CGPL_REPDIAG_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
namespace cgpl
{
	// REPDIAG repeat a matrix along the diagonal a certain number of times, so
	// that if A is a m by n matrix and we want to repeat along the diagonal d
	// times, we get a m*d by n*d matrix B such that:
	// B( (k*m+1):(k*m+1+m-1), (k*n+1):(k*n+1+n-1)) = A 
	// for k from 0 to d-1
	//
	// Inputs:
	//   A  m by n matrix we are repeating along the diagonal. May be dense or
	//     sparse
	//   d  number of times to repeat A along the diagonal
	// Outputs:
	//   B  m*d by n*d matrix with A repeated d times along the diagonal,
	//     will be dense or sparse to match A
	//

	// Sparse version

	template <typename T>
	void repdiag(
		const Eigen::SparseMatrix<T>& A,
		const int d,
		Eigen::SparseMatrix<T>& B);

	template <typename T>
	void repdiag(
		const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A,
		const int d,
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & B);

	// Wrapper with B as output
	template <class Mat>
	Mat repdiag(const Mat & A, const int d);
}

template <typename T>
void cgpl::repdiag(
	const Eigen::SparseMatrix<T>& A,
	const int d,
	Eigen::SparseMatrix<T>& B)
{
	using namespace std;
	using namespace Eigen;
	int m = A.rows();
	int n = A.cols();

	// This will not work for RowMajor
	B.resize(m*d,n*d);
	Eigen::VectorXi per_col = Eigen::VectorXi::Zero(n*d);
	for (int k=0; k<A.outerSize(); ++k)
	{
		for (typename Eigen::SparseMatrix<T>::InnerIterator it(A,k); it; ++it)
		{
			for(int r = 0;r<d;r++) per_col(n*r + k)++;
		}
	}
	B.reserve(per_col);
	for(int r = 0;r<d;r++)
	{
		const int mr = m*r;
		const int nr = n*r;
		for (int k=0; k<A.outerSize(); ++k)
		{
			for (typename Eigen::SparseMatrix<T>::InnerIterator it(A,k); it; ++it)
			{
				B.insert(it.row()+mr,k+nr) = it.value();
			}
		}
	}
	B.makeCompressed();

}

template <typename T>
void cgpl::repdiag(
	const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
	const int d,
	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
	int m = A.rows();
	int n = A.cols();
	B.resize(m*d,n*d);
	B.array() *= 0;
	for(int i = 0;i<d;i++)
	{
		B.block(i*m,i*n,m,n) = A;
	}
}

template <class Mat>
Mat cgpl::repdiag(const Mat & A, const int d)
{
	Mat B;
	repdiag(A,d,B);
	return B;
}
#endif

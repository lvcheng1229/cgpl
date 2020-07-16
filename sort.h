#ifndef CGPL_SORT_H
#define CGPL_SORT_H

#include <vector>
#include <Eigen/Core>
namespace cgpl
{
	template <class T>
	void sort(
		const std::vector<T> & unsorted,
		const bool ascending,
		std::vector<T> & sorted,
		std::vector<size_t> & index_map);

	template <typename DerivedS>
	void sort(
		const Eigen::VectorXd & unsorted,
		const bool ascending,
		Eigen::PlainObjectBase<DerivedS> & sorted,
		Eigen::VectorXi & index_map);

	template<class T> struct IndexLessThan
	{
		IndexLessThan(const T arr) : arr(arr) {}
		bool operator()(const size_t a, const size_t b) const
		{
			return arr[a] < arr[b];
		}
		const T arr;
	};
}

template <typename DerivedS>
void cgpl::sort(
	const Eigen::VectorXd & unsorted,
	const bool ascending,
	Eigen::PlainObjectBase<DerivedS> & sorted,
	Eigen::VectorXi & index_map)
{

	std::vector<double>unsortvec;
	std::vector<double>sortvec;
	std::vector<size_t>indexmap;
	for (int i = 0; i < unsorted.size();i++)
	{
		unsortvec.push_back(unsorted(i));
	}

	sorted.resize(unsorted.rows(), unsorted.cols());
	index_map.resize(unsorted.rows(), unsorted.cols());
	cgpl::sort(unsortvec, ascending, sortvec, indexmap);
	for (int i = 0; i < indexmap.size(); i++)
	{
		sorted(i, 0) = sortvec[i];
		index_map(i, 0) = indexmap[i];
	}
}
template <class T>
void cgpl::sort(
	const std::vector<T> & unsorted,
	const bool ascending,
	std::vector<T> & sorted,
	std::vector<size_t> & index_map)
{
	// Original unsorted index map
	index_map.resize(unsorted.size());
	for(size_t i=0;i<unsorted.size();i++)
	{
		index_map[i] = i;
	}
	// Sort the index map, using unsorted for comparison
	std::sort(
		index_map.begin(),
		index_map.end(),
		cgpl::IndexLessThan<const std::vector<T>& >(unsorted));

	// if not ascending then reverse
	if(!ascending)
	{
		std::reverse(index_map.begin(),index_map.end());
	}
	// make space for output without clobbering
	sorted.resize(unsorted.size());

	// reorder unsorted into sorted using index map
	std::vector<T> copy = unsorted;
	sorted.resize(index_map.size());
	for(int i = 0; i<(int)index_map.size();i++)
	{
		sorted[i] = copy[index_map[i]];
	}
}
#endif
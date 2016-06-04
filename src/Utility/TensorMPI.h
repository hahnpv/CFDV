#pragma once 

	///////////////////////////////////////////////////////////////////////////
	/// TENSOR MPI Parallelized Functions									///
	///////////////////////////////////////////////////////////////////////////

	// Dot product for a 1D tensor (NOT USED IN GMRES)
template<class T>
T MPIdot(Tensor<T, 1> & t1, Tensor<T, 1> & t2)
{
	T result = 0;
	T local_result = dot(t1, t2);

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return result;
}

	// Dot product for a 2D tensor, inner indice
template<class T>
T MPIdot(Tensor<T, 2> & t1, unsigned int j1, Tensor<T, 2> & t2, unsigned int j2)
{
	T result = 0;
	T local_result = dot(t1, j1, t2, j2);

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return result;
}

	// Magnitude of a 1D tensor 
template<class T>
T MPImag(Tensor<T, 1> & t)
{
	T result = 0;
	T local_result = dot(t,t);

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return sqrt(result);
}

	// Magnitude of a 2D tensor, inner column (GMRES)
template<class T>
T MPImag(Tensor<T, 2> & t, int j)
{
	T result = 0;
	T local_result = dot(t, j, t, j);

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return sqrt(result);
}

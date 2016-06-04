#pragma once
//#include "mpi.h"
//#include "add.h"
#include <vector>
#include <ostream>
#include "Utility/split.h"
#include "Utility/tostring.h"

/*
	This version of Tensor uses a std::vector<T> as the storage class.
	It works as-is.
	A template specialization is necessary for T = bool, as bools are 
	stored as bits, not full bytes unless a different memory allocator
	is used. The specialization simply uses an array T* data.
*/

	// imax-lmax are set =1 to multiply to clear data
	// need to NOT DO THAT , instead make a total size variable (size)
	// Can specialize by ndim ...
template< class T, int ndim>
struct Tensor
{
	~Tensor()
	{
	}
		// attempt at copy operator
	template<  class Ti, int ndimi>
	Tensor(const Tensor<Ti, ndimi>& rhs) 
	{
		data = rhs.data;

		imax = rhs.imax;
		jmax = rhs.jmax;
		kmax = rhs.kmax;
		lmax = rhs.lmax;

		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
	}

	Tensor<T, ndim> & operator=(const Tensor<T, ndim>& rhs) 
	{
		data = rhs.data;
	
		imax = rhs.imax;
		jmax = rhs.jmax;
		kmax = rhs.kmax;
		lmax = rhs.lmax;

		x = rhs.x;
		y = rhs.y;
		z = rhs.z;

		return *this;
	}

	Tensor(int lengthi)					// ndim == 1
		: x(0), y(0), z(0),				// i + j*x + k*y + l*z = indice operation if all in the same vec.
		  imax(lengthi), jmax(1), kmax(1), lmax(1)
	{
		data.resize(lengthi);
		for (int i = 0; i < lengthi; i++)
		{
			data[i] = 0;
		}
	}
/*
	Tensor(int lengthi, int lengthj)			// ndim == 2
		: x(lengthi), y(0), z(0),				// i + j*x + k*y + l*z = indice operation if all in the same vec.
 		  imax(lengthi), jmax(lengthj), kmax(1), lmax(1)
	{
		data.resize(lengthi*lengthj);
		for (int i = 0; i < lengthi*lengthj; i++)
		{
			data[i] = 0;
		}
	}
*/
	Tensor(int lengthi, int lengthj, int lengthk)		// ndim == 3
		: x(lengthi), y(lengthi*lengthj), z(0),			// i + j*x + k*y + l*z = indice operation if all in the same vec.
		  imax(lengthi), jmax(lengthj), kmax(lengthk), lmax(1)
	{
		data.resize(lengthi*lengthj*lengthk);
		for (int i = 0; i < lengthi*lengthj*lengthk; i++)
		{
			data[i] = 0;
		}
	}

	Tensor(int lengthi, int lengthj, int lengthk, int lengthl)				// ndim == 4
		: x(lengthi), y(lengthi*lengthj), z(lengthi*lengthj*lengthk),		// i + j*x + k*y + l*z = indice operation if all in the same vec.
		  imax(lengthi), jmax(lengthj), kmax(lengthk), lmax(lengthl)
	{
		data.resize(lengthi*lengthj*lengthk*lengthl);
		for (int i = 0; i < lengthi*lengthj*lengthk*lengthl; i++)
		{
			data[i] = 0;
		}
	}

	T & operator()(int indexi)										// ndim == 1
	{
		return data[indexi];
	}

	T & operator()(int indexi, int indexj)							// ndim == 2
	{
		return data[indexi + x*indexj];
	}
	
	T & operator()(int indexi, int indexj, int indexk)				// ndim == 3
	{
		return data[indexi + x*indexj + y*indexk];
	}

	T & operator()(int indexi, int indexj, int indexk, int indexl)	// ndim == 4
	{
		return data[indexi + x*indexj + y*indexk + z*indexl];
	}

	void clear()
	{
		for (unsigned int i=0; i < imax*jmax*kmax*lmax; i++)
		{
			data[i] = 0;
		}
	}
	unsigned int size( int i)
	{
		if ( i == 0)
			return imax;
		if ( i == 1)
			return jmax;
		if ( i == 2)
			return kmax;
		if ( i == 3)
			return lmax;

//		throw new exception;
	}

	int x, y, z;								// i + j*x + k*y + l*z = indice operation if all in the same vec.
	unsigned int imax, jmax, kmax, lmax;		// think we need a copy operator to make these const.
	
	std::vector<T> data;
};

//////////////////////2D Spec///////////////
template<class T>
struct Tensor <T, 2>
{
	~Tensor()
	{
	}
		// attempt at copy operator
	template<  class Ti, int ndimi>
	Tensor(const Tensor<Ti, ndimi>& rhs) 
	{
		data = rhs.data;

		imax = rhs.imax;
		jmax = rhs.jmax;

		x = rhs.x;
	}

	Tensor<T, 2> & operator=(const Tensor<T, 2>& rhs) 
	{
		data = rhs.data;	

		imax = rhs.imax;
		jmax = rhs.jmax;

		x = rhs.x;

		return *this;
	}

	Tensor(int lengthi, int lengthj)			// ndim == 2
		: x(lengthi), y(0), z(0),				// i + j*x + k*y + l*z = indice operation if all in the same vec.
 		  imax(lengthi), jmax(lengthj), kmax(1), lmax(1)
	{
		data.resize(lengthi*lengthj);
		for (int i = 0; i < lengthi*lengthj; i++)
		{
			data[i] = 0;
		}
	}

	T & operator()(int indexi, int indexj)							// ndim == 2
	{
		return data[indexi + x*indexj];
	}

	void clear()
	{
		for (unsigned int i=0; i < imax*jmax; i++)
		{
			data[i] = 0;
		}
	}
	unsigned int size( int i)
	{
		if ( i == 0)
			return imax;
		if ( i == 1)
			return jmax;

//		throw new exception;
	}

	int x,y,z;								// i + j*x + k*y + l*z = indice operation if all in the same vec.
	unsigned int imax, jmax,kmax,lmax;		// think we need a copy operator to make these const.
	
	std::vector<T> data;
};

////////////////////////////////////////////

template<>
struct Tensor<bool,1>
{
	~Tensor() {};

	Tensor(int lengthi)					// ndim == 1
		: imax(lengthi)
	{
		data = new bool[lengthi];
		for (int i = 0; i < lengthi; i++)
		{
			data[i] = 0;
		}
	}

	 bool & operator()(int indexi)										// ndim == 1
	{
		return data[indexi];
	}

	void clear()
	{
		for (unsigned int i=0; i < imax; i++)
		{
			data[i] = 0;
		}
	}
	unsigned int size( int i)
	{
		if ( i == 0)
			return imax;

//		throw new exception;
	}

	unsigned int imax;	
	
	bool * data;
};

/*
template<> double & Tensor<double, 1>::operator()(int indexi)										// ndim == 1
	{
		return data[indexi];
	}
*/
/*
template<class T> struct Tensor<T, 1>		// attempt at 1D specialization
{
	Tensor(int lengthi)					// ndim == 1
		: x(0), y(0), z(0)				// i + j*x + k*y + l*z = indice operation if all in the same vec.
	{
		imax = lengthi;
		data.resize( lengthi );
	}

	T & operator()(int indexi)										// ndim == 1
	{
		return data[indexi];
	}
};
*/


	///////////////////////////////////////////////////////////////////////////
	/// TENSOR FUNCTIONS													///
	///////////////////////////////////////////////////////////////////////////

	// Magnitude of a 1D tensor
template<class T>
T mag(Tensor<T, 1> & t)
{
	T result = 0;
	for (unsigned int i=0; i < t.imax; i++)
	{
		result += t(i) * t(i);
	}
	return sqrt(result);
}

	// Magnitude of a 2D tensor, inner column
template<class T>
T mag(Tensor<T, 2> & t, int j)
{
	T result = 0;
	for (unsigned int n=0; n<t.imax; n++)
	{	
		result += (t(n, j) * t(n, j));
	}
	return sqrt(result);
}


	// Dot product for a 1D tensor
template<class T>
T dot(Tensor<T, 1> & t1, Tensor<T, 1> & t2)
{
	T result = 0;
	for (int i=0; i < t1.imax; i++)
	{
		result += t1(i) * t2(i);
	}
	return result;
}


	// Dot product for a 2D tensor, inner indice
template<class T>
T dot(Tensor<T, 2> & t1, unsigned int j1, Tensor<T, 2> & t2, unsigned int j2)
{
	T result = 0;
	for (unsigned int i=0; i < t1.imax; i++)
	{
		result += t1(i, j1) * t2(i, j2);
	}
	return result;
}

	// Cross product for 2 1D tensors, size 3.
	// Later, make a generic cross and this a specialized cross.
template<class T>
Tensor<T, 1> cross(Tensor<T, 1> & t1, Tensor<T, 1> & t2)
{
	Tensor<T,1> result(3);
	result(0) = (t1(1)*t2(2) - t1(2)*t2(1));
	result(1) = (t1(2)*t2(0) - t1(0)*t2(2));
	result(2) = (t1(0)*t2(1) - t1(1)*t2(0));
	return result;
}
/*
template<class T>
Tensor<T, 2> Identity(int i)
{
	Tensor<double, 2> Identity(i,i); // temporary?

	Identity.clear();
	for (int ii=0; ii < i; ii++)
	{
		Identity(ii, ii) = 1.0;
	}

	return Identity;
}
*/

template<class T, int ndim>
bool operator==(Tensor<T, ndim> & t1, Tensor<T, ndim> & t2)
{
	for (int i=0; i < t1.imax * t1.jmax * t1.kmax * t1.lmax; i++)
	{
		if ( t1.data[i] != t2.data[i])
			return false;
	}
	return true;
}


// create a function that multiplies all entries together and returns type T.


	///////////////////////////////////////////////////////////////////////////
	/// TENSOR SERIALIZATION												///
	///////////////////////////////////////////////////////////////////////////
	// FORMAT:
	// 1. Type
	// 2. ndim
	// 3. imax, jmax, kmax, lmax. Up to ndim
	// 4. data read straight from data, not in any particular indice order
	// kinda gimp since (1) you need T, ndim to instantiate tensor from the start ...
	// so really you only need the data ... but it works. 
template<class T, int ndim>
std::ostream& write ( std::ostream& os, Tensor<T, ndim> & t )
{
	// lame but all I know to do as of right now
	// or encode in Tensor class... can also typeid().name() but can get fugly
	if (typeid(T) == typeid(double))
	{
		os << "double" << " ";
	}
	else if (typeid(T) == typeid(bool))
	{
		os << "bool" << " ";		
	}
	else if (typeid(T) == typeid(int))
	{
		os << "int" << " ";		
	}

	os << ndim;

	if (ndim > 0)
		os << " " << t.imax;
	if (ndim > 1)
		os << " " << t.jmax;
	if (ndim > 2)
		os << " " << t.kmax;
	if (ndim > 3)
		os << " " << t.lmax;

	for (int i = 0; i < t.imax * t.jmax * t.kmax * t.lmax ; i++)
	{
		os << " " << t.data[i];
	}

	os << endl;

//	cin.get();
	return os;
}

template<class T, int ndim>
std::istream& read ( std::istream& os, Tensor<T, ndim> & t )
{
	std::string line;
	getline(os, line);

	std::vector<std::string> data = split<std::string>(line, " ");

	int count = 0;

	std::string type = data[count++];
/*
	if (type == "double")
	{
		cout << "double" << endl;
	}
	else if (type == "bool")
	{
		cout << "bool" << endl;
	}
	else if (type == "int")
	{
		cout << "int" << endl;
	}
*/
		// should be string_to<int> ?
	int ndimi = string_to<int>(data[count++]);

	int imax = 1, jmax=1, kmax=1, lmax=1;
	if (ndimi > 0)
		imax = string_to<int>(data[count++]);
	if (ndimi > 1)
		jmax = string_to<int>(data[count++]);
	if (ndimi > 2)
		kmax = string_to<int>(data[count++]);
	if (ndimi > 3)
		lmax = string_to<int>(data[count++]);

	for (int i = 0; i < imax * jmax * kmax * lmax ; i++)
	{
		t.data[i] = string_to<T>(data[count++]);
	}

	// FORMAT:
	// 1. Type
	// 2. ndim
	// 3. imax, jmax, kmax, lmax. Up to ndim
	// 4. data read straight from data, not in any particular indice order

	// lame but all I know to do as of right now
	// or encode in Tensor class... can also typeid().name() but can get fugly

//	cin.get();
	return os;
}

	///////////////////////////////////////////////////////////////////////////
	/// TENSOR STREAM PRINTING 												///
	///////////////////////////////////////////////////////////////////////////

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 1> & t ) 
{
	for (int i=0; i < t.imax; i++)
	    os << t(i) << ", ";

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 2> & t ) 
{

	for (int i = 0; i < t.imax; i++)
	{
		for (int j = 0; j < t.jmax; j++)
		{
			os << t(i, j) << ", ";
		}
		os << std::endl;
	}

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 3> & t ) 
{

	for (int i = 0; i < t.imax; i++)
	{
		os << "i = " << i << std::endl;

		for (int j = 0; j < t.jmax; j++)
		{
			for (int k = 0; k < t.kmax; k++)
			{
				os << t(i, j, k) << ", ";
			}
			os << std::endl;
		}
		os << std::endl;
	}

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 4> & t ) 
{

	for (int i = 0; i < t.imax; i++)
	{
		for (int j = 0; j < t.jmax; j++)
		{
			os << "(i, j) = " << i << ", " << j << std::endl;

			for (int k = 0; k < t.kmax; k++)
			{
				for (int l = 0; l < t.lmax; l++)
				{
					os << t(i, j, k) << ", ";
				}
				os << std::endl;
			}
			os << std::endl;
		}
	}

	return os;
}

	///////////////////////////////////////////////////////////////////////////
	/// TENSOR MPI Parallelized Functions									///
	///////////////////////////////////////////////////////////////////////////


	// Dot product for a 1D tensor (NOT USED IN GMRES)
template<class T>
T MPIdot(Tensor<T, 1> & t1, Tensor<T, 1> & t2)
{
	T result = 0;
	T local_result = 0;
	for (int i=0; i < t1.imax; i++)
	{
		local_result += t1(i) * t2(i);
	}

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return result;
}

	// Dot product for a 2D tensor, inner indice (GMRES)
template<class T>
T MPIdot(Tensor<T, 2> & t1, unsigned int j1, Tensor<T, 2> & t2, unsigned int j2)
{
	T result = 0;
	T local_result = 0;
	for (unsigned int i=0; i < t1.imax; i++)
	{
		local_result += t1(i, j1) * t2(i, j2);
	}

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return result;
}

	// Magnitude of a 1D tensor (GMRES)
	// appears to work (both threads return same value)
	// might try mpi_reduce_scatter?
template<class T>
T MPImag(Tensor<T, 1> & t)
{
	T result = 0;
	T local_result = 0;

//	cout << MPI::COMM_WORLD.Get_rank() << " summing ... " << endl;
	for (unsigned int i=0; i < t.imax; i++)
	{
		local_result += t(i) * t(i);
	}
//	cout << MPI::COMM_WORLD.Get_rank() << " sum complete, reducing ... " << local_result << endl;

//	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

//	cout << MPI::COMM_WORLD.Get_rank() << " Allreduce complete " << endl;
	return sqrt(result);
}

	// Magnitude of a 2D tensor, inner column (GMRES)
template<class T>
T MPImag(Tensor<T, 2> & t, int j)
{
	T result = 0;
	T local_result = 0;
	for (unsigned int n=0; n<t.imax; n++)
	{	
		local_result += (t(n, j) * t(n, j));
	}

	MPI::COMM_WORLD.Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM);

	return sqrt(result);
}

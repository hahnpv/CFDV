#pragma once

#include <cmath>
#include <vector>
#include <boost/shared_ptr.hpp>
// using namespace boost;

template<class T>
struct RefCnt
{
	RefCnt(int lengthi, int lengthj, int lengthk, int lengthl)
		: x(lengthi),
		  y(lengthi*lengthj),
		  z(lengthi*lengthj*lengthk),
		  imax(lengthi), 
		  jmax(lengthj), 
		  kmax(lengthk), 
		  lmax(lengthl)
	{
		size = imax + x*jmax + y*kmax + z*lmax;
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = 0;
	}

	~RefCnt()
	{
		delete [] data;
	}

	void clear()
	{
		for (int i=0; i < size; i++)
		{
			data[i] = 0;
		}
	}


	T * data;
	int size;
	int x, y, z;								// strides
	unsigned int imax, jmax, kmax, lmax;		// think we need a copy operator to make these const.
};


//enum indice { i, j, k, l, m, n, o, p, q, r, s, t };

	// Can and should specialize by ndim ... so ndim==2 cant use operator(i,j,k), etc.
	// and should re-write MPI operations in terms of the standard operations, etc.
	// Separate functionality into files like MPI, ostream, etc.

template< class T, int ndim>
struct Tensor
{
	Tensor(const Tensor<T, ndim> & rhs)									// CopyOp
	{
		ref = rhs.ref;
	}

	Tensor<T, ndim> & operator=(const Tensor<T, ndim> & rhs)			// AssignmentOp
	{
		ref = rhs.ref;
		return *this;
	}

	Tensor(int lengthi)													// ndim == 1
		: ref( new RefCnt<T>(lengthi, 0, 0, 0) )
	{}

	Tensor(int lengthi, int lengthj)									// ndim == 2
		: ref( new RefCnt<T>(lengthi, lengthj, 0, 0) )
	{}

	Tensor(int lengthi, int lengthj, int lengthk)						// ndim == 3
		: ref( new RefCnt<T>(lengthi, lengthj, lengthk, 0) )
	{}

	Tensor(int lengthi, int lengthj, int lengthk, int lengthl)			// ndim == 4
		: ref( new RefCnt<T>(lengthi, lengthj, lengthk, lengthl) )
	{}

		  // At some point in time, specialize by ndim.
	T & operator()(int indexi)										// ndim == 1
	{
		return ref->data[indexi];
	}

	T & operator()(int indexi, int indexj)							// ndim == 2
	{
		return ref->data[indexi + ref->x*indexj];
	}
	
	T & operator()(int indexi, int indexj, int indexk)				// ndim == 3
	{
		return ref->data[indexi + ref->x*indexj + ref->y*indexk];
	}

	T & operator()(int indexi, int indexj, int indexk, int indexl)	// ndim == 4
	{
		return ref->data[indexi + ref->x*indexj + ref->y*indexk + ref->z*indexl];
	}

	void clear()
	{
		ref->clear();
	}

	unsigned int imax()	{ return ref->imax; };
	unsigned int jmax()	{ return ref->jmax; };
	unsigned int kmax()	{ return ref->kmax; };
	unsigned int lmax()	{ return ref->lmax; };

//	std::vector< indice > index;						// indice

	boost::shared_ptr< RefCnt<T> > ref;
};

/*
// implement for only 1dim tensors ... can't?
template<class T, int ndim> T & Tensor<T, ndim>::operator()(int indexi)						// ndim == 1
{
	return ref->data[indexi];
}

template<class T, int ndim> T & Tensor<T, 1>::operator()(int indexi)						// ndim == 1
{
	return ref->data[indexi];
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

	// Dot product for a 1D tensor
template<class T>
T dot(Tensor<T, 1> & t1, Tensor<T, 1> & t2)
{
	T result = 0;
	for (unsigned int i=0; i < t1.imax(); i++)
	{
		result += t1(i) * t2(i);
	}
	return result;
}

	// Magnitude of a 1D tensor
template<class T>
T mag(Tensor<T, 1> & t)
{
	return std::sqrt( dot( t, t) );
}

	// Dot product for a 2D tensor, inner indice
template<class T>
T dot(Tensor<T, 2> & t1, unsigned int j1, Tensor<T, 2> & t2, unsigned int j2)
{
	T result = 0;
	for (unsigned int i=0; i < t1.imax(); i++)
	{
		result += t1(i, j1) * t2(i, j2);
	}
	return result;
}

	// Magnitude of a 2D tensor, inner column
template<class T>
T mag(Tensor<T, 2> & t, int j)
{
	return std::sqrt( dot(t,j, t,j) );
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

template<class T,int i>
Tensor<T, 2> Kroneker_delta()
{
	return Identity<T>(i);
}


template<class T, int ndim>
bool operator==(Tensor<T, ndim> & t1, Tensor<T, ndim> & t2)
{
	for (int i=0; i < t1.imax() * t1.jmax() * t1.kmax() * t1.lmax(); i++)
	{
		if ( t1.ref->data[i] != t2.ref->data[i])
			return false;
	}
	return true;
}

	// sum all elements in tensor
template<class T, int ndim>
T sum(Tensor<T,ndim> & t)
{
	int sum = 0;
	for (int i = 0; i < t.ref->size(); i++)
	{
		sum += t.ref->data[i];
	}
	return sum;
}
	// aggregate (mult) all members
	// better way would be stl algorithm route, eventually.
template<class T, int ndim>
T aggregate(Tensor<T,ndim> & t)
{
	T sum = 1.0;
	for (int i = 0; i < t.ref->size; i++)
	{
		sum *= t.ref->data[i];
	}

	return sum;
}

Tensor<double, 1> unit(Tensor<double,1> & t)
{
	Tensor<double, 1> retval( t.imax() );
	double length = mag(t);
	for (unsigned int i = 0; i < t.imax(); i++)
		retval(i) = t(i) / length;
	return retval;
}


// create a function that multiplies all entries together and returns type T.

#include "TensorMath.h"					// Mathematics
#include "TensorMPI.h"					// MPI parallelized functionale
#include "TensorStream.h"				// std::ostream printing

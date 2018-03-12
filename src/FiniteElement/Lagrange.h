#pragma once

#include "../Utility/Tensor.h"

	// Lagrange interpolation Polynomial (nondimensionalized [0,1])
	// USAGE:
	//	Lagrange L(basis);
	// for (int i = 0; i < 10; i++)
	//     result += L[0.55](i) * line(i);
struct Lagrange : public Tensor<double, 1>
{
	Lagrange(Tensor<double, 1> & base)		// number of samples
		: Tensor<double,1>(base.imax()),
		  basis(base.imax()),
		  size(base.imax()-1)
	{
		for (int i = 0; i <= size; i++)			
		{
			basis(i) = base(i);				// copy data
		}
	}

	Lagrange & operator[](double xi)		// uniform distribution
	{
		if (xi == xi_0)				// don't re-calc; allows L[N](i) usage.
			return *this;

		xi_0 = xi;

		for (int N = 0; N <= size; N++)
		{
			(*this)(N)	= 1;

			for (int M = 0; M <= size; M++)
			{
				if ( M != N)
				{
					(*this)(N) *= ( xi - basis(M) ) / ( basis(N) - basis(M)  );
				}
			}
		}
		return *this;
	}

	double xi_0;
	Tensor<double, 1> basis;
	unsigned int size;
};

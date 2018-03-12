#pragma once

// tau functional object
#include <iostream>
#include "../Utility/Tensor.h"
using namespace std;
struct Tau : public Tensor<double, 2>
{
	Tau(int dimn) : Tensor<double, 2>(dimn, dimn) // should automatically scale from 2D to 3D. Validate each case.
	{
		dim = dimn;
	}

	void calculate(double mu, Tensor<double, 2> & dV)
	{
		clear();

		Tensor<double, 2> delta = Kroneker_delta<double,2>();

		double vkk = 0;

		for (unsigned int k=0; k < dV.imax(); k++)
		{
			vkk += dV(k, k);
		}

		for (int i=0; i < dim; i++)
		{
			for (int j=0; j < dim; j++)
			{
				(*this)(i, j) = mu * (dV(i, j) + dV(j, i));

				if ( i == j)
				{
					(*this)(i, j) -= (2.0/3.0) * mu * vkk;				// -lambda * vkk
				}
		
				// could be re-expressed as
		//		(*this)(i, j) = mu * (dV(i, j) + dV(j, i) - (2.0/3.0) * vkk * delta(i, j));
			}
		}
	}

	int dim;
};

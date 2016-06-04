#pragma once

#include "tau.h"

// flow is a container for holding variables in domain and boundary integral
// so I can instantiate one flow in element, and have all the memory allocated permenantly
// since same vars are used in domain and boundary, and no persistance, can use one struct.

struct Flow
{
	Flow(int ndim, int neqn)
		: 	v(ndim),
			x(ndim),
			U(neqn),
			dU(neqn,ndim),
			dV(ndim,ndim),
			dT(ndim),
			tau(ndim)
	{}

	void clear()
	{
		E  = 0;
		T  = 0;
		v.clear();
		x.clear();
		U.clear();
		dU.clear();
		dV.clear();
		dT.clear();
		k  = 0;
		mu = 0;
		tau.clear();
	}

	void dump()
	{
	}

	double E;
	double T;
	Tensor<double, 1> v;
	Tensor<double, 1> x;
	Tensor<double, 1> U;
	Tensor<double, 2> dU;
	Tensor<double, 2> dV;
	Tensor<double, 1> dT;
	double k;
	double mu;
	Tau tau;
};


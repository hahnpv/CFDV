#pragma once

#include "Utility/Tensor.h"

	// struct for the values at each node in an element
	// think about implementing inheritance, a base node class
	// and then bc node class(es): either a single one or multiple
	// then can filter out (have a vector of pointers to bc nodes in main node vector)
	// the bc node classes and just send them for dirichlet bc application
	// small optimization but in a million node mesh might make a difference
struct Node
{
	Node(int ndim, int neqn)
		: v(ndim),
		  U(neqn),
		  U0(neqn),
		  x(ndim),
		  dirichlet(neqn)
	{
		for (unsigned int i=0; i<dirichlet.imax(); i++)
			dirichlet(i) = false;
	}

		// Flow Variables
	double rho;									/// Density
	Tensor<double, 1> v;						/// Velocity vector
	double E;									/// Energy
	double e;									/// E - 0.5 * V^2

	double p;									/// Pressure
	double T;									/// Temperature

		// Turbulence Variables (k-w)
	double k;
	double w;

		// Navier-Stokes conservation variable
	Tensor<double, 1> U;						/// Navier-Stokes Conservation Variable
	Tensor<double, 1> U0;						/// Old value of the conservation variable
	
		// Geometry
	Tensor<double, 1> x;						/// Location
	unsigned int number;						/// node number

		// Boundary Conditions
//	int bc;
	Tensor<bool, 1> dirichlet;

	int c;
	double coeff;
};


#pragma once

#include "FDV/Node.h"

// 1D test functions for line integrals (boundaries)

struct Phibou : public Tensor<double, 1> // Matrix1d
{
	Phibou(Tensor<double, 1> & zeta) : Tensor<double, 1>(2) // Matrix1d(2)
	{	
		(*this)(0) = 0.5 * (1 - zeta(0));
		(*this)(1) = 0.5 * (1 + zeta(0));
	}
};

struct dPxi : public Tensor<double, 1> // Matrix1d
{
	dPxi() : Tensor<double, 1>(2) // Matrix1d(2)
	{
		(*this)(0) = -0.5;
		(*this)(1) =  0.5;
	}
};

struct Normal : public Tensor<double, 1> // Matrix1d
{
	Normal(Tensor<double, 1> & w, /*dPxi dpxi,*/ Tensor<double, 1> & dpxi, std::vector<Node *>node) : Tensor<double, 1>(2) // Matrix1d(2)
	{
		double dxds = 0;
		double dyds = 0;
		for (int i = 0; i < 2; i++)						// iterate over both nodes
		{
			dxds += dpxi(i) * node[i]->x(0);
			dyds += dpxi(i) * node[i]->x(1);
		}
		(*this)(0) =  dyds * w(0);
		(*this)(1) = -dxds * w(0);
	}
};

struct dPidx : Tensor<double, 2> // Matrix2d
{
	dPidx(Node & left, Node & right) : Tensor<double, 2>(2,2) // Matrix2d(2,2)
	{
		double delx = right.x(0) - left.x(0);
		double dely = right.x(1) - left.x(1);

		double dpdx[2] = {0, 0};
		double dpdy[2] = {0, 0};
		if (abs(delx) < 1.0E-12)
		{}							// =0
		else
		{
			dpdx[0] = -1.0/delx;
			dpdx[1] =  1.0/delx;
		}
		if (abs(dely) < 1.0E-12)
		{}							// =0
		else
		{
			dpdy[0] = -1.0/dely;
			dpdy[1] =  1.0/dely;
		}

		for (int i = 0; i < 2; i++)			// iterate over both nodes
		{
			(*this)(i, 0) = dpdx[i];
			(*this)(i, 1) = dpdy[i];
		}
	}
};

struct TFNorm2D //: public Matrix1d
{
	TFNorm2D(std::vector<Node *> & node, Tensor<double, 1> & phibou, Tensor<double, 1> & en, Tensor<double, 2> & dshpdx, Tensor<double, 1> & w, Tensor<double, 1> & zeta) //: Matrix1d(3)
	{
		phibou = Phibou(zeta);
		dPxi dpxi;
		en = Normal(w,dpxi,node);
		dshpdx = dPidx(*node[0],*node[1]);
	}
};

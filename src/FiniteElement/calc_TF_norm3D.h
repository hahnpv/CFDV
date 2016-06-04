#pragma once

#include "FiniteElement/TestFunctions2D.h"					// dpxi

	// calculate the test function and normal to a 3D boundary surface necessary for integration

struct TFNorm3D 
{
	TFNorm3D(std::vector<Node *> & node, Tensor<double, 1> & phibou, Tensor<double, 1> & en, Tensor<double, 2> & dshpdx, Tensor<double, 1> & w, Tensor<double, 1> & zeta, Element * e)
	{
		int nnod = node.size();			// number of nodes, 4 on the boundary
		int ndim = 3;					// 3D

		// phibou = testfunction(zeta);
		phibou(0) = 0.25 * (1-zeta(0)) * (1-zeta(1));
		phibou(1) = 0.25 * (1+zeta(0)) * (1-zeta(1));
		phibou(2) = 0.25 * (1+zeta(0)) * (1+zeta(1));
		phibou(3) = 0.25 * (1-zeta(0)) * (1+zeta(1));

		Tensor<double, 1> Lambda(ndim);			// vector e12
		Tensor<double, 1> Mu(ndim);				// vector e14

		Tensor<double, 2> nodeloc(4,3);
		nodeloc.clear();
		if (!e->adap)
		{
			for (int N = 0; N < 4; N++)
			{
				for (int j = 0; j < 3; j++)
				{
					nodeloc(N, j) = node[N]->x(j);
				}
			}
		} 
		else
		{
			Tensor<double, 2> allnode(8,3);
			allnode.clear();
		
			for (int N = 0; N < e->Hnm.imax(); N++)
			{
				for (int M = 0; M < e->Hnm.jmax(); M++)
				{
					for (int j = 0; j < 3; j++)
					{
//						dPxiHnm(M, j) += dPxi(N, j) * Hnm(N, M);
						allnode(N, j) += e->Hnm(N, M) * e->node[M]->x(j);
					}
				}
			}

			// now find nodes in allnode to put in nodeloc
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (node[j]->number == e->node[i]->number)
					{
						for (int k = 0; k < 3; k++)
						{
							nodeloc(j, k) = allnode(i, k);
						}
					}
				}
			}
		}

		for(int i = 0; i < ndim; i++)
		{
			Lambda(i) = nodeloc(1,    i) - nodeloc(0,    i);
			Mu(i)     = nodeloc(3,    i) - nodeloc(0,    i);
		}

		for(int j = 0; j < ndim; j++)
		{
			Lambda(j) /= mag(Lambda);							// each mag is calculated 3 times, should probably cache result
			Mu(j)     /= mag(Mu);
		}

//		Tensor<double, 1> enorm = cross(Lambda, Mu);		// vector normal to surface
//		Tensor<double, 1> gama = cross(enorm, Lambda);		// y-prime axis
		Tensor<double, 1> enorm(3);
		Tensor<double, 1> gama(3);
	enorm(0) = (Lambda(1)*Mu(2) - Lambda(2)*Mu(1));
	enorm(1) = (Lambda(2)*Mu(0) - Lambda(0)*Mu(2));
	enorm(2) = (Lambda(0)*Mu(1) - Lambda(1)*Mu(0));

	gama(0) = (enorm(1)*Lambda(2) - enorm(2)*Lambda(1));
	gama(1) = (enorm(2)*Lambda(0) - enorm(0)*Lambda(2));
	gama(2) = (enorm(0)*Lambda(1) - enorm(1)*Lambda(0));

		// coordinate transformation matrix
		Tensor<double, 2> acoortr(ndim,ndim);
		for( int i = 0; i < ndim; i++)
		{
			acoortr(0, i) = Lambda(i);
			acoortr(1, i) = gama(i);
			acoortr(2, i) = enorm(i);
		}

		// components of surface Jacobian
		dtestfunction dpxi(zeta);

		// could transform phi, dphi here ? then use node[N]->x(i)
		// in theory, identical, test out later.

		Tensor<double, 2> dx_dxi(3,2);							// surface derivatives
		dx_dxi.clear();
        for( int N = 0; N < 4; N++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					dx_dxi(i, j) += dpxi(N, j) * /*node[N]->x(i);*/ nodeloc(N, i);
				}
			}
		 }

			// Transformed 2D Jacobian
		Tensor<double, 2> dxdxitr(2,2);							// surface derivatives transformed to xyz coord space
		dxdxitr.clear();
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					dxdxitr(i, j) += acoortr(i, k) * dx_dxi(k, j);
				}
			}
		}

		double detjac = dxdxitr(0, 0) * dxdxitr(1, 1) - dxdxitr(1, 0) * dxdxitr(0, 1);
        double djacinv = 1./detjac;

		if ( detjac < 1.0E-12)
		{
			cout << "Face Normal surface integral has a negative detjac: " << detjac << endl;
			cout << nodeloc << endl;
			cout << acoortr << endl;
			cout << "in element " << e->number << endl;
			cin.get();
		}

		// Calculate the derivative with respect to zetas
		Tensor<double, 2> dpdx(nnod,2);
        for( int N = 0; N < 4; N++)
		{
			dpdx(N, 0) = ( dxdxitr(1, 1) * dpxi(N, 0) - dxdxitr(1, 0) * dpxi(N, 1)) * djacinv;
		    dpdx(N, 1) = (-dxdxitr(0, 1) * dpxi(N, 0) + dxdxitr(0, 0) * dpxi(N, 1)) * djacinv;
		}

		// Calculate the spatial derivatice of the shape function and the normal vector
		for (int i=0; i < ndim; i++)
		{
			for (int N = 0; N < nnod; N++)
			{
				for (int j = 0; j < 2; j++)
				{
					dshpdx(N, i) += acoortr(j, i) * dpdx(N, j);							// dpidx
				}
			}

            en(i) = enorm(i) * w(0) * w(1) * detjac;									// normal
		}
	}
	
};


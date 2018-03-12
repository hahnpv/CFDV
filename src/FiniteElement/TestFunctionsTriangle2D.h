#pragma once

#include "../Utility/Tensor.h"
#include "../FDV/Node.h"

	// test function, 2d, triangle, per Chung ch. 9.
struct testfunctiontri : public Tensor<double, 1> // Matrix1d
{
	testfunctiontri(Tensor<double, 1> & zeta) : Tensor<double, 1>(3) // Matrix1d(4)
	{
		(*this)(0) = 0.25 * (1-zeta(0)) * (1-zeta(1));
		(*this)(1) = 0.25 * (1+zeta(0)) * (1-zeta(1));
		(*this)(2) = 0.5 * (1+zeta(1));
	}
};

	// derivative of the test function WRT natural coordinates
struct dtestfunctiontri : public Tensor<double, 2> // Matrix2d
{
	dtestfunctiontri(Tensor<double, 1> & zeta) : Tensor<double, 2>(3, 2) // Matrix2d(4,2)
	{
		(*this)(0, 0) = -0.25 * (1-zeta(1));		// d/xi
		(*this)(1, 0) =  0.25 * (1-zeta(1));
		(*this)(2, 0) =  0;

		(*this)(0, 1) = -0.25 * (1-zeta(0));		// d/deta
		(*this)(1, 1) = -0.25 * (1+zeta(0));
		(*this)(2, 1) =  0.5;
	}
};

struct jacobiantri	: Tensor<double, 2> // Matrix2d		// eventaully put dx/dy dxi/deta into an array can automate more calculations, etc.
{
	jacobiantri(Tensor<double, 2> & dPxi, std::vector<Node *> & node) : Tensor<double, 2>(2, 2) // Matrix2d(2,2)
	{
		clear();
		for (int N = 0; N < 3; N++)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					(*this)(i, j) += dPxi(N, j) * node[N]->x(i);
				}
			}
		}

		J = ((*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1));
	}

	double det()
	{
		return J;
	};

	double J;
};

struct dtestfunctiondxtri : public Tensor<double, 2> // Matrix2d
{
	dtestfunctiondxtri(Tensor<double, 2> & dPxi, jacobiantri & j) : Tensor<double, 2>(3, 2) // Matrix2d(4,2)
	{
		double dxdxi  = j(0, 0);
		double dxdeta = j(0, 1);
		double dydeta = j(1, 1);
		double dydxi  = j(1, 0);

		double J = j.det();

		for(int N=0; N < 3; N++)
		{
			(*this)(N, 0) = (1.0/J) * (dPxi(N, 0)*dydeta - dydxi*dPxi(N, 1));
			(*this)(N, 1) = (1.0/J) * (dxdxi*dPxi(N, 1) - dPxi(N, 0)*dxdeta);
		}		
	}
};
  
/*
struct d2testfunction : public Tensor<double, 1> // Matrix1d
{
	d2testfunction() : Tensor<double, 1>(3) // Matrix1d(4)
	{
		(*this)(0) =  0.25;
		(*this)(1) = -0.25;
		(*this)(2) =  0.25;		// ?
	}
};
*/
/*
struct d2phidx2 : public Matrix3d
{
	d2phidx2( std::vector<Node *> & node, jacobian j, dtestfunction dpxi, d2testfunction d2phi ) : Matrix3d(2,2,4)
	{
		// extract Jacobian data
		double dydeta = j.dydeta;
		double dydxi  = j.dydxi;
		double dxdxi  = j.dxdxi;
		double dxdeta = j.dxdeta;
		double detjac = j.J;
		// extract Jacobian data

		double dxxieta = 0;
		double dyxieta = 0;
		for( int inode=0; inode<4; inode++)
		{
			dxxieta += d2phi[inode] * node[inode]->x[0];
			dyxieta += d2phi[inode] * node[inode]->x[1];
		}
		double djacdxi = dxdxi*dyxieta - dxxieta*dydeta;
		double djacdeta = dxxieta*dydeta - dxdeta*dyxieta;
		double detjac2 = 1./pow(detjac,2);								// fortran evaluates exponents before mult/div so should just modify J
		double detjac3 = 1./pow(detjac,3);

		double term1 = 0;
		double term2 = 0;
		double term3 = 0;
		double term4 = 0;
//		Matrix3d d2phidx(2,2,4);
		for( int inode=0; inode<4; inode++)
		{
			term1 = -detjac3*djacdxi*(pow(dydeta,2)*dpxi[inode, 0] - dydeta*dydxi*dpxi[inode, 1]);	
			term2 = detjac2*(dydeta*dyxieta*dpxi[inode, 0] - dydeta*dydxi*d2phi[inode]);
			term3 = detjac3*djacdeta*(dydxi*dydeta*dpxi[inode, 0] - pow(dydxi,2)*dpxi[inode, 1]);
			term4 = -detjac2*(dydxi*dydeta*d2phi[inode] - dydxi*dyxieta*dpxi[inode, 1]);


		//!     d2phi/dx2
			(*this)(0, 0, inode) = term1 + term2 + term3 + term4;

			term1 = detjac3*djacdxi*(-pow(dxdeta,2)*dpxi[inode, 0] + dxdeta*dxdxi*dpxi[inode, 1]);
			term2 = -detjac2*(-dxdeta*dxxieta*dpxi[inode, 0]  + dxdeta*dxdxi*d2phi[inode]);
			term3 = -detjac3*djacdeta*(-dxdxi*dxdeta*dpxi[inode, 0] + pow(dxdxi,2)*dpxi[inode, 1]);
			term4 = detjac2*(-dxdxi*dxdeta*d2phi[inode] + dxdxi*dxxieta*dpxi[inode, 1]);

		//!     d2phi/dy2
			(*this)(1, 1, inode) = term1 + term2 + term3 + term4;

			term1 = -detjac3*djacdxi*(-dydeta*dxdeta*dpxi[inode, 0] + dydeta*dxdxi*  dpxi[inode, 1]);
			term2 = detjac2*(-dydeta*dxxieta*dpxi[inode, 0]+ dydeta*dxdxi*d2phi[inode]);
			term3 = detjac3*djacdeta*(-dydxi*dxdeta*dpxi[inode, 0] + dydxi*dxdxi*dpxi[inode, 1]);
			term4 = -detjac2*(-dydxi*dxdeta*d2phi[inode] + dydxi*dxxieta*dpxi[inode, 1]);

		//!     d2phi/dxdy
			(*this)(0, 1, inode) = term1 + term2 + term3 + term4;

			term1 = detjac3*djacdxi*(dxdeta*dydeta*dpxi[inode, 0] - dxdeta*dydxi*dpxi[inode, 1]);
			term2 = -detjac2*(dxdeta*dyxieta*dpxi[inode, 0] - dxdeta*dydxi*d2phi[inode]);
			term3 = -detjac3*djacdeta*(dxdxi*dydeta*dpxi[inode, 0] - dxdxi*dydxi*dpxi[inode, 1]);
			term4 = detjac2*(dxdxi*dydeta*d2phi[inode] - dxdxi*dyxieta*dpxi[inode, 1]);

		//!     d2phidydx
			(*this)(1, 0, inode) = term1 + term2 + term3 + term4;
		}
	}
};
*/

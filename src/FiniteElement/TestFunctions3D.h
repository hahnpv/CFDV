#pragma once

#include "FDV/Node.h"
// #include "MeshRefine2D.h"

	// test function, 3D
struct testfunction3d : Tensor<double, 1> // Matrix1d
{
//	~testfunction3d() { std::cout << " ~testfunction3d " << std::endl; }
	testfunction3d(Tensor<double, 1> & zeta) : Tensor<double, 1>(8) // Matrix1d(8)
	{
		(*this)(0) = 0.125 * (1-zeta(0)) * (1-zeta(1)) * (1-zeta(2));
		(*this)(1) = 0.125 * (1+zeta(0)) * (1-zeta(1)) * (1-zeta(2));
		(*this)(2) = 0.125 * (1+zeta(0)) * (1+zeta(1)) * (1-zeta(2));
		(*this)(3) = 0.125 * (1-zeta(0)) * (1+zeta(1)) * (1-zeta(2));
		(*this)(4) = 0.125 * (1-zeta(0)) * (1-zeta(1)) * (1+zeta(2));
		(*this)(5) = 0.125 * (1+zeta(0)) * (1-zeta(1)) * (1+zeta(2));
		(*this)(6) = 0.125 * (1+zeta(0)) * (1+zeta(1)) * (1+zeta(2));
		(*this)(7) = 0.125 * (1-zeta(0)) * (1+zeta(1)) * (1+zeta(2));
	}
};

	// derivative of the test function WRT natural coordinates
struct dtestfunction3d : public Tensor<double, 2> // Matrix2d
{
//	~dtestfunction3d() { std::cout << " ~dtestfunction3d " << std::endl; }
	dtestfunction3d(Tensor<double, 1> & zeta) : Tensor<double, 2>(8, 3) // Matrix2d(8,3)
	{
		(*this)(0, 0) = -0.125 * (1-zeta(1)) * (1-zeta(2));		// d/xi
		(*this)(1, 0) =  0.125 * (1-zeta(1)) * (1-zeta(2));
		(*this)(2, 0) =  0.125 * (1+zeta(1)) * (1-zeta(2));
		(*this)(3, 0) = -0.125 * (1+zeta(1)) * (1-zeta(2));
		(*this)(4, 0) = -0.125 * (1-zeta(1)) * (1+zeta(2));
		(*this)(5, 0) =  0.125 * (1-zeta(1)) * (1+zeta(2));
		(*this)(6, 0) =  0.125 * (1+zeta(1)) * (1+zeta(2));
		(*this)(7, 0) = -0.125 * (1+zeta(1)) * (1+zeta(2));

		(*this)(0, 1) = -0.125 * (1-zeta(0)) * (1-zeta(2));		// d/deta
		(*this)(1, 1) = -0.125 * (1+zeta(0)) * (1-zeta(2));
		(*this)(2, 1) =  0.125 * (1+zeta(0)) * (1-zeta(2));
		(*this)(3, 1) =  0.125 * (1-zeta(0)) * (1-zeta(2));
		(*this)(4, 1) = -0.125 * (1-zeta(0)) * (1+zeta(2));
		(*this)(5, 1) = -0.125 * (1+zeta(0)) * (1+zeta(2));
		(*this)(6, 1) =  0.125 * (1+zeta(0)) * (1+zeta(2));
		(*this)(7, 1) =  0.125 * (1-zeta(0)) * (1+zeta(2));

		(*this)(0, 2) = -0.125 * (1-zeta(0)) * (1-zeta(1));		// d/dzeta
		(*this)(1, 2) = -0.125 * (1+zeta(0)) * (1-zeta(1));
		(*this)(2, 2) = -0.125 * (1+zeta(0)) * (1+zeta(1));
		(*this)(3, 2) = -0.125 * (1-zeta(0)) * (1+zeta(1));
		(*this)(4, 2) =  0.125 * (1-zeta(0)) * (1-zeta(1));
		(*this)(5, 2) =  0.125 * (1+zeta(0)) * (1-zeta(1));
		(*this)(6, 2) =  0.125 * (1+zeta(0)) * (1+zeta(1));
		(*this)(7, 2) =  0.125 * (1-zeta(0)) * (1+zeta(1));
	}
};
	// eventually matrixify it, gettign convoluted. Probably universal from 2D to 3D
struct jacobian3d : public Tensor<double, 2> // Matrix2d			// eventaully put dx/dy dxi/deta into an array can automate more calculations, etc.
{
//	~jacobian3d() { std::cout << " ~jacobian3d " << std::endl; }
	jacobian3d(Tensor<double, 2> & dPxi, std::vector<Node *> & node) : Tensor<double, 2>(3, 3) // Matrix2d(3,3)
	{
		clear();
		for (unsigned int N=0; N<node.size(); N++)		// removed j=0, clear(), should be redundant.
		{
			for (int i = 0; i < 3; i++)						// dx, dy, dx
			{
				for (int j = 0; j < 3; j++)					// dxi, deta, dzeta
				{
					(*this)(i, j) += dPxi(N, j) * node[N]->x(i);
				}
			}
		}

		J = (*this)(0, 0)*(*this)(1, 1)*(*this)(2, 2) 
		  + (*this)(0, 1)*(*this)(1, 2)*(*this)(2, 0) 
		  + (*this)(0, 2)*(*this)(1, 0)*(*this)(2, 1) 
		  - (*this)(0, 2)*(*this)(1, 1)*(*this)(2, 0) 
		  - (*this)(0, 1)*(*this)(1, 0)*(*this)(2, 2) 
		  - (*this)(0, 0)*(*this)(1, 2)*(*this)(2, 1);
	}

	double det()
	{
		return J;
	};

	double J;
};

struct dtestfunctiondx3d : public Tensor<double, 2> //  Matrix2d
{
//	~dtestfunctiondx3d() { std::cout << " ~dtestfunction3ddx " << std::endl; }
	dtestfunctiondx3d(Tensor<double, 2> & dPxi, jacobian3d & j) : Tensor<double, 2>(8, 3) // Matrix2d(8,3)
	{
		// invert jacobian
//		Matrix2d jinv(3,3);
		Tensor<double, 2> jinv(3, 3);

		jinv(0, 0)=   j(1, 1)*j(2, 2) - j(1, 2)*j(2, 1);
		jinv(0, 1)= -(j(1, 0)*j(2, 2) - j(1, 2)*j(2, 0));
		jinv(0, 2)=   j(1, 0)*j(2, 1) - j(1, 1)*j(2, 0);
		jinv(1, 0)= -(j(0, 1)*j(2, 2) - j(0, 2)*j(2, 1));
		jinv(1, 1)=   j(0, 0)*j(2, 2) - j(0, 2)*j(2, 0);
		jinv(1, 2)= -(j(0, 0)*j(2, 1) - j(0, 1)*j(2, 0));
		jinv(2, 0)=   j(0, 1)*j(1, 2) - j(0, 2)*j(1, 1);
		jinv(2, 1)= -(j(0, 0)*j(1, 2) - j(0, 2)*j(1, 0));
		jinv(2, 2)=   j(0, 0)*j(1, 1) - j(0, 1)*j(1, 0);
		
		double J = j.det();

		for (int i=0; i < 3; i++)
		{
			for(int N=0; N < 8; N++)
			{
				for(int j=0; j < 3; j++)
				{
					(*this)(N, i) += (1.0/J) * ( jinv(i, j) * dPxi(N, j) );
				}
			}
		}
 	}
};

	// TEMPORARY until TF's are unified or some other course of action taken.
  void tf3d(Tensor<double, 1> & Phi, double & J, Tensor<double, 2> & dPhi, Tensor<double, 1> & zeta, std::vector<Node *> & node)
  {
  		Phi = testfunction3d(zeta);
		dtestfunction3d dPxi(zeta);
		jacobian3d jac(dPxi, node);
		dPhi = dtestfunctiondx3d(dPxi,jac);

		J = jac.det();
  }

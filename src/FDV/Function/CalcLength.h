#pragma once

#include <iostream>
using namespace std;

#include "../../FiniteElement/get_face.h"
#include "../../FiniteElement/TestFunctionsTriangle2D.h"

	// Quadrature and Test Functions
#include "../../FiniteElement/GaussQuad.h"
#include "../../FiniteElement/TestFunctions2D.h"
#include "../../FiniteElement/TestFunctions3D.h"

	// Adaptive Meshing
#include "../../adap/MeshRefine2D.h"			/// adaptive

	/// Calculate the characteristic length of each element
template<class T> struct CalcLength : public unary_function<T, void>
{
	CalcLength(int ndim, int nnod)
		: gq(2),
		  ndim(ndim),
		  nnod(nnod)
	{}

	~CalcLength()
	{
		e = 0;			// move pointer from a legitimate element, or Bad Things Happen
		delete e;
		e = NULL;
	}

	typedef void result_type;					// result_type for tr1 wrapper.
	result_type operator() (T e) 
	{
		this->e = e;

		e->area = 0;

//		ele_t eletype = get_ele_t(ndim,nnod);
		if (ndim == 2 && nnod == 3)
		{
			gq.two(this,&CalcLength::updateTri);
		}
		else if (ndim == 2)
		{
			gq.two(this,&CalcLength::update2D);
		}
		else if (ndim == 3)
		{
			gq.three(this,&CalcLength::update3D);
		}
		else if (ndim == 1)
		{
			e->area = abs( e->node[1]->x(0) - e->node[0]->x(0) );
		}

		e->cleng = pow(e->area,1.0/(double)ndim);


		if ( e->area < 0)
		{	
			cout << "area: " << e->area << " " << e->cleng << " in FE " << e->number << endl;
			cin.get();
		}
	}

		/// 2D triangle finite element area integration
	void updateTri(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
		dtestfunctiontri dPxi(zeta);

		jacobiantri jac(dPxi, e->node);
		double J = jac.det();

		e->area += w(0) * w(1) * J;
	}

		/// 2D finite element area integration
	void update2D(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
		dtestfunction dPxi(zeta);
	
		Tensor<double, 1> dummy(4);
		if (e->adap)
		{
		//	Adaptive::adap( Tensor<double, 1>(4), dPxi, e->Hnm);
		//	adap( dummy, dPxi, e->Hnm);
		}

		jacobian jac(dPxi, e->node);
		double J = jac.det();
		e->area += w(0) * w(1) * J;
	}

		/// 3D finite element area integration
	void update3D(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
		dtestfunction3d dPxi(zeta);
		Tensor<double, 1> dummy(8);

		if (e->adap)
		{
//			adap( dummy, dPxi, e->Hnm);
		}

		jacobian3d jac(dPxi, e->node);
		double J = jac.det();
		e->area += w(0) * w(1) * w(2) * J;
	}

	int ndim;
	int nnod;
	Element * e;
	GaussQuad<CalcLength> gq;
};


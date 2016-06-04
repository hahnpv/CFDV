#pragma once

#include <vector> 

#include "FiniteElement/Face.h"
#include "Node.h"
#include "Utility/Tensor.h"

enum ele_t { iso2d, iso3d, tria2d, iso1d };		// element type, isoparametric or triangular

	// the struture for an element - should be idential in header
struct Element
{
	Element(const int nnod, const int neqn, const int ndim)			// need to verify dimensioning on ALL mats
	:	R(nnod*neqn,nnod*neqn),
		rhs(nnod*neqn),
		L(nnod*neqn,nnod*neqn),
		U(nnod*neqn,nnod*neqn),
		Hnm(nnod, nnod),
		adap(false),
		refine_level(0)
{
	// Nodal values
	node.resize( nnod);
}
	int number;							/// element number
	double cleng;						/// characteristic length
	double area;						/// area of the finite element

	std::vector< Node *> node;			/// Flow variables at each node (F, G)
	std::vector<Face *> face;			/// Boundary Conditions	(was <Face>)

	Tensor<double,2> R;					/// Left hand side of the system of equations (A + B)
	Tensor<double,1> rhs;				/// Right hand side of the system of equations (H + N)

	Tensor<double, 2> L;				/// Lower matrix, LU Factorization
	Tensor<double, 2> U;				/// Upper matrix, LU Factorization

	Tensor<double, 2> Hnm;				/// Adaptive meshing translation matrix
	bool adap;							/// adaptive flag, true if Hnm is active.
	int refine_level;					/// refine level of the element in question.
	ele_t type;							/// element type (not in use yet)
};

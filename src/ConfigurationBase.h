#pragma once

// Base class for CFD configurations
// for instance,
// -> unparallelized
// -> MPI
// -> CUDA
// etc...
#include "Utility/dictionary.h"
#include "FDV/Element.h"
#include "Utility/Solvers/GMRESbase.h"
struct ConfigurationBase
{
	ConfigurationBase(dictionary & key)
	{
		nnod = key.get_val<int>("nnod");							/// Number of nodes per finite element
		neqn = key.get_val<int>("neqn");							/// Number of equations per node
		ndim = key.get_val<int>("ndim");							/// Number of dimensions 
	}

	~ConfigurationBase()
	{
		delete gmres;
	}

	virtual void reset_gmres() = 0;

	virtual void Output(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, double time, bool ele_data) = 0;

	virtual void RMSErr(double t, int iter) = 0;

	virtual void Save(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, dictionary & key) = 0;
	
	int rank;

	int nnod;
	int neqn;
	int ndim;
	std::vector<Node *>::iterator NodeIteratorStart;
	std::vector<Node *>::iterator NodeIteratorEnd;
	GMRESbase * gmres;
};
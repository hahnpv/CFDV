#pragma once

#include "../../FDV/Element.h"
#include "../../FDV/Thermo.h"
#include "../split.h"
#include "../tostring.h"
#include "../dictionary.h"

#include <iostream>
#include <fstream>

	// Base class for all data loaders, defines common functionality
struct Load
{
	Load(std::vector<std::string> & args,
		 std::vector< Element * > & elements, 
		 std::vector< Node *> & nodes
		 ): elements(elements),
			nodes(nodes)
	{
	}

protected:


		/// Reads a generic tagged input file
		/// key = val
	void read_file(std::string filename)
	{
		ifstream file(filename.c_str(),ios::in);
		if (!file.is_open())
		{
			cout << "Error opening " << filename << "!" << endl;
			cin.get();
			return;
		}
		else
		{
			while(!file.eof())
			{
				string input;
				getline(file, input);
				key.add_key(split<string>(input," \t"));
			}
		}
		file.close();
	}

		/// Add a node to the node list
	void add_node(std::vector<double>  node)
	{
		Thermo & thermo = Thermo::Instance();
	
		if (node.size() < 2*ndim+3)
			return;

		int i = (int)node[0] - 1;		// human count, not c++ count, fix consistency

		nodes[i] = new Node(ndim,neqn);

		for (int j=0; j < ndim; j++)
		{
			nodes[i]->x(j)	= node[j+1];
		}

		// Build the conservation variable
		for(int j=0; j<ndim+2; j++)
		{
			nodes[i]->U(j) = node[(unsigned int)ndim+1+j];
		}

		// conversion from E input to T input:
		nodes[i]->U((unsigned int)ndim+1) *= thermo.Cv;					// *= Cv

		double v = 0;
		for (int j=0; j<ndim; j++)
		{
			v += nodes[i]->U(j+1) * nodes[i]->U(j+1);					// += 0.5*V^2i		// could use dot outside loop
			nodes[i]->U(j+1) *= nodes[i]->U(0);							// sum rho into velocity
		}

		nodes[i]->U((unsigned int)ndim+1) += 0.5 * v;					// sum 1/2*V^2
		nodes[i]->U((unsigned int)ndim+1) *= nodes[i]->U(0);			// times rho

//		double bc = node[2*(unsigned int)ndim+3];		// number,x,y,rho,u,v,T, bc=7 for 2D
														// number,x,y,z,rho,u,v,w,T, bc=9 for 3D -> 2*ndim+3
//		nodes[i]->bc   = (int)bc;						// currently only time its used is to flat Adiabatic in applybc

		nodes[i]->number = i;											// deprecate this.
	}

		/// Add an element to the element list
	void add_element(std::vector<int>  elem)
	{
		if (elem.size() < nnod+1)
			return;

		int i = elem[0] - 1;			// human count, not c++ count, need to fix consistency.

		elements[i] = new Element((int)nnod, (int)neqn, (int)ndim);
		elements[i]->number = i;

		for (int j=0; j < nnod; j++)
		{
			elements[i]->node[j] = nodes[elem[j+1]-1];
		}
	}

		/// Add a face to an element
	void add_face(std::vector<int> face)	// take a face line and apply to elems
	{
		if ( face.size() == 0)
			return;

	//	cout << "face.size(): " << face.size() << endl;

		Thermo & thermo = Thermo::Instance();

	//	cout << "nbnod = " << nbnod << endl; 

		Face * f = new Face(nbnod);				// should be nnod/2, number of nodes

		// Need to deprecate the face variable slot in  the input file ... its null.

		f->face = face[1];
		f->bc   = face[2];

		for (int i=0; i < nbnod; i++)
		{
			f->n(i) = (unsigned int)face[i+3];
		}

		int bc = f->bc;

		elements[ face[0] ]->face.push_back( f);		// push a pointer instead of make a copy

			cout << "el: " << face[0] << " nodes: " << f->n << endl;
		for (int i=0; i< nbnod; i++)
		{
			if ( bc == -1 || bc == -12 )		// inflow, lid driven flow
			{
				for( unsigned int j=0; j<nodes[i]->dirichlet.imax(); j++)
				{
					elements[ face[0] ]->node[f->n(i)]->dirichlet(j) = true;
				}
			}

			else if ( bc == 2 || bc == 9)	// reflection
			{
				elements[ face[0] ]->node[f->n(i)]->dirichlet(ndim) = true;
			}

			// TEST - not sure ... 
			// reflection in XZ plane, mirroring reflection in the XY plane, not sure, check book.
			else if ( bc == 22)
			{
				elements[ face[0] ]->node[f->n(i)]->dirichlet(2) = true;	
			}

			else if ( bc == 1 || bc == 11 || bc == 21 || bc == 31 || bc == 41 || bc == 4 )	// wall, compresion corner
			{
				// 1     - XY flat plate
				// 11,21 - YZ flat plate
				// 31,41 - XZ flat plate
				// 4     - ramp
				for( int j=0; j<ndim; j++)
				{
					elements[ face[0] ]->node[f->n(i)]->dirichlet(1+j) = true;
				}
				if( !thermo.adiabatic)
				{
					elements[ face[0] ]->node[f->n(i)]->dirichlet(1+ndim) = true;
				}
			}
		}
	}

		/// Create the thermodynamic property structure
	void add_thermo()
	{
		Thermo & thermo = Thermo::Instance();

		thermo.setGamma( key.get_val<double>("gamma"));	
		thermo.Pr    = key.get_val<double>("Pr");
		thermo.csuth = 110.0 / key.get_val<double>("Tinf");
		thermo.cmach = key.get_val<double>("M");
		thermo.Cv    = 1./thermo.gamma/(thermo.gamma-1.0)/pow(thermo.cmach, 2);
		thermo.cgas  = 1.0/1.4/pow(thermo.cmach, 2);
		thermo.creyn = key.get_val<double>("Re");
		thermo.Twall = key.get_val<double>("Twall");
		thermo.adiabatic = key.get_val<bool>("Adiabatic");
	}

	std::vector< Element * > & elements;
	std::vector< Node *> & nodes;

	std::string aerofile;
	std::string contfile;
	std::string gridfile;

	int nnod;
	int neqn;
	int ndim;
	int nbnod;

public:
	dictionary key;						// make this a public dictionary for the CFD 
};

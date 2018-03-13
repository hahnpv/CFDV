#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "../../../FDV/Thermo.h"
#include "../../../FiniteElement/get_face.h"
#include "../../dictionary.h"
#include "../../split.h"

#include "../../../FDV/Element.h"

#include "boost/program_options.hpp"
namespace po = boost::program_options;


using namespace std;

	/// 3D
struct node_3d_t
{
	double x[3];
	double rho;
	double v[3];
	double T;
};

struct ele_3d_t
{
	int node[8];
	int bc[6];
};

	/// 2D
struct node_2d_t
{
	double x[2];
	double rho;
	double v[2];
	double T;
};

struct ele_2d_t
{
	int node[4];
	int bc[4];
};

	/// 1D
struct node_1d_t
{
	double x[1];
	double rho;
	double v[1];
	double T;
};
struct ele_1d_t 
{
	int node[2];
	int bc[2];
};


ostream & operator<<(ostream & os, ele_3d_t &ele)
{
	os << "nodes: ";
	for (int i = 0; i < 8; i++)
		cout << ele.node[i] << " ";
	os << "  bc: ";
	for (int i = 0; i < 6; i++)
		cout << ele.bc[i] << " ";
	
	return os;
}
ostream & operator<<(ostream & os, ele_2d_t &ele)
{
	os << "nodes: ";
	for (int i = 0; i < 4; i++)
		cout << ele.node[i] << " ";
	os << "  bc: ";
	for (int i = 0; i < 4; i++)
		cout << ele.bc[i] << " ";
	
	return os;
}
ostream & operator<<(ostream & os, node_3d_t &n)
{
	os << "coord: " << n.x[0] << " " << n.x[1] << " " << n.x[2] << " cvar: " << n.rho << " " << n.v[0] << " " << n.v[1] << " " << n.v[2] << " " << n.T;
	return os;
}
ostream & operator<<(ostream & os, node_2d_t &n)
{
	os << "coord: " << n.x[0] << " " << n.x[1] << " cvar: " << n.rho << " " << n.v[0] << " " << n.v[1] << " " << n.T;
	return os;
}

struct LoadBinaryData
{
	std::string path;
		
	LoadBinaryData(
			 std::vector< Element * > & elements,
			 std::vector< Node *> & nodes,
	  	     po::variables_map am,
			 std::string path,
			 int nnod, int neqn, int ndim, int nbnod
		)
			 :	elements( elements),
				nodes( nodes),
		        path(path),
				nnod(nnod),
				neqn(neqn),
				ndim(ndim),
				nbnod(nbnod)
	{
		add_thermo(am);
	}
	
	void read_elements(int min, int max, std::string fname = "//binaryelements")
	{
		if ( ndim == 3)
		{
			loadElements<ele_3d_t>(path + fname, min, max);
		}
		else if (ndim == 2)
		{
			loadElements<ele_2d_t>(path + fname, min, max);		
		}
		else if (ndim == 1)
		{
			loadElements<ele_1d_t>(path + fname, min, max);		
		}
	}

	int nodemin;
	void read_nodes(int min, int max, std::string fname = "//binarynodes")
	{
		nodemin = min;	// fixme i dont like this

		if ( ndim == 3)
		{
			loadNodes<node_3d_t>(path + fname, min, max);
		}
		else if (ndim == 2)
		{
			loadNodes<node_2d_t>(path + fname, min, max);
		}
		else if (ndim == 1)
		{
			loadNodes<node_1d_t>(path + fname, min, max);
		}
	}

	void read_adap()
	{
		adap(path + "//adap");
	}

	template< class node_t>
	void loadNodes(std::string path, int min, int max)
	{
		Thermo & thermo = Thermo::Instance();

		ifstream nodein(path.c_str(),ios::binary);

		int offset = sizeof(node_t) * min;

		nodein.seekg( offset, ios::cur);

		for (int i = min; i <= max; i++)
		{
			node_t node;
			nodein.read((char*)&node, sizeof(node_t));

			Node * n = new Node(ndim,neqn);

			for (int j=0; j < ndim; j++)
			{
				n->x(j)	= node.x[j];
			}

				/// Build the conservation variable, U
			n->U(0) = node.rho;
			for (int j = 0; j < ndim; j++)
			{
				n->U(j+1) = node.v[j];
			}

			n->U(ndim+1) = node.T * thermo.Cv;

			double v = 0;
			for (int j=0; j<ndim; j++)
			{
				v += n->U(j+1) * n->U(j+1);					// += 0.5*V^2i		// could use dot outside loop
				n->U(j+1) *= n->U(0);							// sum rho into velocity
			}

			n->U((unsigned int)ndim+1) += 0.5 * v;					// sum 1/2*V^2
			n->U((unsigned int)ndim+1) *= n->U(0);			// times rho

			n->number = i;											// deprecate this.

			nodes.push_back( n);
		}

		nodein.close();
	};

	template< class ele_t>
	void loadElements(std::string path, int min, int max)
	{
		ifstream elein(path.c_str(),ios::binary);
		int nface = 6;	// 3D
		if (ndim == 2)
			nface = 4;
		else if (ndim == 1)
			nface = 2;

		int offset = sizeof(ele_t) * min;
		elein.seekg( offset, ios::cur);

		for (int i = min; i < max; i++)
		{
			ele_t ele;

			elein.read((char*)&ele, sizeof(ele_t));

			Element * e = new Element((int)nnod, (int)neqn, (int)ndim);
			e->number = i;

			for (int j=0; j < nnod; j++)
			{
				e->node[j] = nodes[ele.node[j] - nodemin]; // nodes[elem[j+1]-1];
			}

			for (int j = 0; j < nface; j++)
			{
				if (ele.bc[j] != 0)
				{
					add_face( e, j, ele.bc[j]);
				}
			}

			elements.push_back( e);
		}

		elein.close();
	};

	void adap(std::string path)
	{
		std::string input;
		ifstream file(path.c_str());
		while (!file.eof())
		{
			getline(file,input);
//			add_adap( split<double>(input," \t") );
			std::vector<double> data = split<double>(input, " \t");
			add_adap( data);
		}
	};

	void add_adap(std::vector<double> & data)
	{
		// 0 - element
		// 1 - adap flag
		// 2 - tensor data
		if (data.size() == 0)
			return;

		int elno = (int)data[0];
		bool adap = (bool)data[1];
		int refine_level = (int)data[2];

		elements[ elno ]->adap = adap;
		elements[ elno ]->refine_level = refine_level;

		std::string datastring = "double ";
		for (unsigned int i = 4; i < data.size(); i++)			// this all is hackish, get working then fix.
		{												// might make more sense to make a string ver. of read
			datastring += to_string<double>(data[i]);	// and have the stream version call it once it stringsplits
			datastring += " ";
		}
		std::istringstream is(datastring);

		read( is, elements[ elno ]->Hnm);
	}

		/// Add a face to an element
	void add_face(Element * e, int face, int bc)	// take a face line and apply to elems
	{
		Thermo & thermo = Thermo::Instance();

		Face * f = new Face(nbnod);				// should be nnod/2, number of nodes

		f->face = face;
		f->bc   = bc;

		f->n = get_nodes_from_face( face, get_ele_t(ndim,nnod) );

		e->face.push_back( f);		// push a pointer instead of make a copy

		for (int i=0; i< nbnod; i++)
		{
			if ( bc == -1 || bc == -12 )		// inflow, lid driven flow
			{
				for( unsigned int j=0; j<nodes[i]->dirichlet.imax(); j++)
				{
					e->node[f->n(i)]->dirichlet(j) = true;
				}
			}

			else if ( bc == 2 || bc == 9)	// reflection
			{
				e->node[f->n(i)]->dirichlet(ndim) = true;		// rho = 0, u = 1, v = 2, w = 3
			}

			// TEST - not sure ... 
			// reflection in XZ plane, mirroring reflection in the XY plane, not sure, check book.
			else if ( bc == 22)
			{
				e->node[f->n(i)]->dirichlet(2) = true;	
			}

			else if ( bc == 1 || bc == 11 || bc == 21 || bc == 31 || bc == 41 || bc == 4 )	// wall, compresion corner
			{
				// 1     - XY flat plate
				// 11,21 - YZ flat plate
				// 31,41 - XZ flat plate
				// 4     - ramp
				cout << "i = " << i << endl;
				cout << "BC: node " << f->n(i) << endl;
				cout << "dirichlet size: " << e->node[f->n(i)]->dirichlet.imax() << endl;
				for( int j=0; j<ndim; j++)
				{
					e->node[f->n(i)]->dirichlet(1+j) = true;
				}
				if( !thermo.adiabatic)
				{
					e->node[f->n(i)]->dirichlet(1+ndim) = true;
				}
			}
		}
	}
		/// Create the thermodynamic property structure
	void add_thermo(po::variables_map am)
	{
		Thermo & thermo = Thermo::Instance();

		thermo.setGamma( am["gamma"].as<double>());
		thermo.Pr    = am["Pr"].as<double>();
		thermo.csuth = 110.0 / am["Tinf"].as<double>();
		thermo.cmach = am["M"].as<double>();
		thermo.Cv    = 1./thermo.gamma/(thermo.gamma-1.0)/pow(thermo.cmach, 2);
		thermo.cgas  = 1.0/1.4/pow(thermo.cmach, 2);
		thermo.creyn = am["Re"].as<double>();
		thermo.Twall = am["Twall"].as<double>();
		thermo.adiabatic = am["Adiabatic"].as<bool>();
	}

	std::vector< Element * > & elements;
	std::vector< Node *> & nodes;

	int nnod;
	int neqn;
	int ndim;
	int nbnod;
};

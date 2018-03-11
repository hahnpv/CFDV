#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

#include "../../../FDV/Thermo.h"
#include "../../../FiniteElement/get_face.h"
#include "../../dictionary.h"
#include "../../split.h"

#include "../../../FDV/Element.h"

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
			 std::vector<std::string> & args,
			 std::vector< Element * > & elements, 
			 std::vector< Node *> & nodes 
			 )
			 :	elements( elements),
				nodes( nodes)
		{
	
		// args points to DIRECTORY instead of FILES.
		if (args.size() < 2)			// catch debug instantiation
		{
		//	path = "..//input//flat_plate_validate_fortran_new";
			path = "input//3dtestlagrange";
		//	path = "input//smallbackstep";
		}
		else
		{
			path = args[1];
		}

		key.add_key("path", path);

		cout << "path: " << path << endl;
		cout << "aeropath: " << path + "//aero" << endl;

		aero(path + "//aero");
		control(path + "//control");

		// set up nnod, neqn, ndim
		nnod = key.get_val<int>("nnod");
		neqn = key.get_val<int>("neqn");
		ndim = key.get_val<int>("ndim");
		nbnod = key.get_val<int>("nbnod");

		// need to set up thermo before nodes, for Cv in E calc.
		add_thermo();	
	}
	
	void read_elements(int min, int max, std::string fname = "//binaryelements")
	{
			// set up node/element size here eventually -> in MPI from criteron.
	//	elements.resize( /*key.get_val<int>("elements")*/ max - min );

		// Load nodes, elements, faces. Add parameters for min and max. 
		// Call MPIConfiguration prior to this and use parameters for our rank.
		// 
		if ( ndim == 3)
		{
			cout << "3D" << endl;
//			loadNodes<node_3d_t>(path + "//nodes");
			loadElements<ele_3d_t>(path + fname, min, max);
		}
		else if (ndim == 2)
		{
			cout << "2D" << endl;
//			loadNodes<node_2d_t>(path + "//nodes");
			loadElements<ele_2d_t>(path + fname, min, max);		
		}
		else if (ndim == 1)
		{
			cout << "1D" << endl;
//			loadNodes<node_1d_t>(path + fname, min, max);
			loadElements<ele_1d_t>(path + fname, min, max);		
		}
	}

	int nodemin;
	void read_nodes(int min, int max, std::string fname = "//binarynodes")
	{
		nodemin = min;
			// set up node/element size here eventually -> in MPI from criteron.
	//	nodes.resize( /*key.get_val<int>("nodes")*/ max - min );

		// Load nodes, elements, faces. Add parameters for min and max. 
		// Call MPIConfiguration prior to this and use parameters for our rank.
		// 
		if ( ndim == 3)
		{
			cout << "3D" << endl;
			loadNodes<node_3d_t>(path + fname, min, max);
//			loadElements<ele_2d_t>(path + "//elements");
		}
		else if (ndim == 2)
		{
			cout << "2D" << endl;
			loadNodes<node_2d_t>(path + fname, min, max);
//			loadElements<ele_2d_t>(path + "//elements");		
		}
		else if (ndim == 1)
		{
			cout << "1D" << endl;
			loadNodes<node_1d_t>(path + fname, min, max);
//			loadElements<ele_1d_t>(path + "//elements");		
		}
	}

	void read_adap()
	{
		if ( key.get_val<bool>("adap") )
		{
			cout << "adap" << endl;
			adap(path + "//adap");
		}
	}

	void aero(std::string path)
	{
		read_file(path);
	};

	void control(std::string path)
	{
		read_file(path);
	};

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
			/*
			// HAX missing BC 
			if (i == 1409)
			{
				cout << "HAX IN LoadBinary.h" << endl;
				ele.bc[3] = -1.;
			}
			*/
//			cout << "ele.node: " << ele.node[0] << " " << ele.node[1] << endl;
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

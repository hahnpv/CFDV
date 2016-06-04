#pragma once

#include <iostream>
#include <fstream>

using namespace std;

#include "LoadBinary.h"
#include "FiniteElement/get_face.h"


template< class ele_t, class node_t>
struct SaveBinaryData
{
	SaveBinaryData(
			 std::vector< Element * > & elements, 
			 std::vector< Node *> & nodes,
			 int ndim,
			 int neqn,
			 int nnod
			 )
			 : ndim(ndim),
			   neqn(neqn),
			   nnod(nnod),
			   elements(elements),
			   nodes(nodes)
	{
		cout << "ndim = " << ndim << endl;

		nface = 6;			// 3D
		if ( ndim == 2)
			nface = 4;		// 2D
		if ( ndim == 1)
			nface = 2;		// 1D
	}

	void write_elements(std::string fname)
	{
		ofstream elout(fname.c_str(),ios::binary);
		for (int i = 0; i < elements.size(); i++)
		{
			ele_t ele;
			cout << "nnod/nface/face size: " << nnod << " " << nface << " " << elements[i]->face.size() << endl;
			for (int j = 0; j < nnod; j++)
			{
				ele.node[j] = elements[i]->node[j]->number;
			}
			for (int j = 0; j < nface; j++)
			{
				ele.bc[j] = 0;
			}
			for (int j = 0; j < elements[i]->face.size(); j++)
			{
				int faceno = get_face( elements[i]->face[j]->n, get_ele_t(ndim,nnod) );

				cout << "face = " << faceno << " bc: " << elements[i]->face[j]->bc << " n: " << elements[i]->face[j]->n << endl;
				
				int bc = elements[i]->face[j]->bc;
				ele.bc[ faceno ] = bc;

			}
			cout << "nodes: " << ele.node[0] << " " << ele.node[1] << " bc: " << ele.bc[0] << " " << ele.bc[1] << endl;

			elout.write((char*)&ele, sizeof(ele_t));
		}
		elout.close();
	}

	void write_nodes(std::string fname)
	{
		ofstream nout(fname.c_str(),ios::binary);
		for (int i = 0; i < nodes.size(); i++)
		{
			node_t node;
			for (int j = 0; j < ndim; j++)
			{
				node.x[j] = nodes[i]->x(j);
			}
			node.rho = nodes[i]->rho;
			for (int j = 0; j < ndim; j++)
			{
				node.v[j] = nodes[i]->v(j);
			}
			node.T = nodes[i]->T;

	//		cout << i << " " << node << endl;
			nout.write((char*)&node, sizeof(node_t));
		}
		nout.close();
	}

	int ndim;
	int neqn;
	int nnod;
	int nface;
	std::vector< Element * > elements;
	std::vector< Node *> nodes;
};
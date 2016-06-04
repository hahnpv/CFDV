#pragma once

#include <iostream>
using namespace std;

	/// Class for playing around with mesh refinement concepts
	/// needs to be re-worked completely, start on paper and also
	/// revisit Gary's dissertation and code.
	/// clean and functionalize, make easy to understand.
	/// conceptually it's relatively simple, splitting 1 element into 4 pieces
	/// try to convey that in the code.
	/// may need some kind of a list object to hang on to transformation history.
template<class T> struct MeshRefine : public unary_function<T, void>
{
	MeshRefine(std::vector<T> & elements, std::vector<Node *> & nodes,
				int ndim, int neqn, int nnod)		// dumb, find a better way
		: gq(2),
		  zeta(5,2),	
		  elements(elements),
		  nodes(nodes),
		  ndim(ndim),
		  neqn(neqn),
		  nnod(nnod),
		  phi(nnod)
	{
		zeta[0][0] =  0.0;		// face 0
		zeta[0][1] = -1.0;

		zeta[1][0] =  1.0;		// face 1
		zeta[1][1] =  0.0;

		zeta[2][0] =  0.0;		// face 2
		zeta[2][1] =  1.0;

		zeta[3][0] = -1.0;		// face 4
		zeta[3][1] =  0.0;

		zeta[4][0] =  0.0;		// middle
		zeta[4][1] =  0.0;
	}

	Matrix1d phi;

		// Formulate test function for the point of interest
	Matrix1d & interp(int point)
	{
		phi.clear();

		phi[0] = 0.25 * (1-zeta[point][0]) * (1-zeta[point][1]);
		phi[1] = 0.25 * (1+zeta[point][0]) * (1-zeta[point][1]);
		phi[2] = 0.25 * (1+zeta[point][0]) * (1+zeta[point][1]);
		phi[3] = 0.25 * (1-zeta[point][0]) * (1+zeta[point][1]);

		return phi;
	}

	void operator() (T e) 
	{
		std::vector< int> node;						// temporary storage, global node numbers of local nodes

		for (int i=0; i<4; i++)
		{
			node.push_back( e->node[i]->number);
		}

		Thermo & thermo = Thermo::Instance();

		if(e->s1 > 0.99)
		{
			// need to create 5 nodes
			for (int i=0; i<5; i++)
			{
				phi = interp(i);		// redundant

				cout << "creating node ... " << endl;
				nodes.push_back( new NodeVars(thermo.gamma,thermo.Cv,ndim,neqn));

				int ii = nodes.size() - 1;

					// interpolate x,y
				for (int N=0; N<e->node.size(); N++)
				{
					nodes[ii]->x[0] += phi[N] * e->node[N]->x[0];
					nodes[ii]->x[1] += phi[N] * e->node[N]->x[1]; 


						// interpolate conservation variable
					for (int j=0; j<4; j++)
					{
						nodes[ii]->U[j] += phi[N] * e->node[N]->U[j];
					}
				}

				nodes[ii]->number = ii;

				node.push_back( ii );

				// Determine node BC
				if ( i == 0)
				{
//					if (e->node[0]->bc == e->node[1]->bc)
//					{
						nodes[ii]->bc = e->node[0]->bc;			// along bottom should always be safe to follow left, confirm.
						for (int j = 0; j < e->node[0]->dirichlet.size(); j++)
						{
							nodes[ii]->dirichlet[j] = e->node[0]->dirichlet[j];
						}
//					}
				}
				else if ( i == 1)
				{
					if (e->node[1]->bc != 0 || e->node[2]->bc != 0)
					{
						if (e->node[1]->bc == e->node[2]->bc)
						{
							nodes[ii]->bc = e->node[1]->bc;
							for (int j = 0; j < e->node[1]->dirichlet.size(); j++)
							{
								nodes[ii]->dirichlet[j] = e->node[1]->dirichlet[j];
							}
						}
						else
						{
							cout << "vertical boundary, does not match" << e->node[1]->bc << " " << e->node[1]->bc << endl;
							nodes[ii]->bc = 0;
						}
					}
					else
					{
						nodes[ii]->bc = 0;
					}

				}
				else if ( i == 2)
				{
					nodes[ii]->bc = e->node[0]->bc;
				}
				else if ( i == 3)
				{
					if (e->node[3]->bc != 0 || e->node[0]->bc != 0)
					{
						if (e->node[3]->bc == e->node[0]->bc)
						{
							nodes[ii]->bc = e->node[3]->bc;
							for (int j = 0; j < e->node[3]->dirichlet.size(); j++)
							{
								nodes[ii]->dirichlet[j] = e->node[3]->dirichlet[j];
							}
						}
						else
						{
							cout << "vertical boundary, does not match" << e->node[3]->bc << " " << e->node[0]->bc << endl;
						}
					}
					else
					{
						nodes[ii]->bc = 0;
					}
				}
				else if ( i == 4)
				{
					nodes[ii]->bc = 0;
				}
			}

			// create three additional elements (4 total)
			// define elements 2,3,4
			// (nodes, numbers, faces)
			for (int i=0; i<3; i++)
			{
				cout << "creating element ... " << endl;
				elements.push_back( new element(nnod, neqn, ndim));

				int ii = elements.size() - 1;

				elements[ii]->number = ii;

				if ( i == 0)
				{
					elements[ii]->node[0] = nodes[ node[4] ];				
					elements[ii]->node[1] = nodes[ node[1] ];
					elements[ii]->node[2] = nodes[ node[5] ];
					elements[ii]->node[3] = nodes[ node[8] ];

					for (int j=0; j<e->face.size(); j++)
					{
						if (e->face[j].face == 0 || e->face[j].face == 1)
						{
							elements[ii]->face.push_back( e->face[j] );
						}
					}
				}
				else if ( i == 1)
				{
					elements[ii]->node[0] = nodes[ node[8] ];
					elements[ii]->node[1] = nodes[ node[5] ];
					elements[ii]->node[2] = nodes[ node[2] ];
					elements[ii]->node[3] = nodes[ node[6] ];

					for (int j=0; j<e->face.size(); j++)
					{
						if (e->face[j].face == 1 || e->face[j].face == 2)
						{
							elements[ii]->face.push_back( e->face[j] );
						}
					}
				}
				else if ( i == 2)
				{
					elements[ii]->node[0] = nodes[ node[7] ];				
					elements[ii]->node[1] = nodes[ node[8] ];
					elements[ii]->node[2] = nodes[ node[6] ];
					elements[ii]->node[3] = nodes[ node[3] ];
					
					for (int j=0; j<e->face.size(); j++)
					{
						if (e->face[j].face == 2 || e->face[j].face == 3)
						{
							elements[ii]->face.push_back( e->face[j] );
						}
					}
				}
			}

			// re-define element 1
//			e->node[0] = nodes[ node[0] ];		// same				
			e->node[1] = nodes[ node[4] ];
			e->node[2] = nodes[ node[8] ];
			e->node[3] = nodes[ node[7] ];
			for (int j=0; j<e->face.size(); j++)
			{
				if (e->face[j].face == 1 || e->face[j].face == 2)
				{
					vector<Face>::iterator iter = e->face.begin() + j;
					e->face.erase(iter);		// figure out safe way
				}
			}

			// regenerate area, cleng after mods.
			// (for now on whole domain just cause its easy, later optimize and just do updated elements)
		}

		// unrefine -> e->s1 < 0.2, but need to look at neighbors
	}

	void update(Matrix1d & w, Matrix1d & zeta)
	{

	}

	Element * e;					// should base off template
	GaussQuad<MeshRefine> gq;		// 2D
	Matrix2d zeta;
	std::vector<T> & elements;
	std::vector<Node *> & nodes;

	int ndim, neqn, nnod;

};
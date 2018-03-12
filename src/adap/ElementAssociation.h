#pragma once

// this is the new EleCrunch


	// try using inheritance? face boundary / face elemnts
template<class Element>
struct face
{
	face()
	{
		boundary = false;
	}
	std::vector<Element *> e;	// element(s) on the face
	bool boundary;				// is a domain boundary, MUTUALLY EXCLUSIVE to face.
//	if (!boundary, elements.size())
};

template<class Element>
struct ele
{
	Element * e;							// the element we are
	std::vector< face<Element> > f;			// faces of the element
	std::vector< face<Element> > edge;
};


template<class Element, class Node>
struct node
{
	Node * n;					// the node we are
	std::vector<Element *> e;	// elements at the node
	// global or local correlation?
};

	// This works great but gets slow as things get big.
	// We can probably intelligently kill things early when faces/elements fill up to capacities.
	// also watch when/how often it's called. Pass around instead of re-instantiating.
template<class Elements, class Nodes>
struct ElementAssociation
{
	ElementAssociation(std::vector<Elements *> & e, std::vector<Nodes *> & n, int nnod, int nface)
		: nnod(nnod),
		  nface(nface)
	{
			// call reset to perform initialization 
			// one interface to reset when needed.
		reset(e, n);
	};

	void reset(std::vector<Elements *> & e, std::vector<Nodes *> & n)
	{
		elements.clear();
		nodes.clear();

		elements.resize( e.size() );
		nodes.resize( n.size() );

		cout << "node correlation" << endl;

		// 0. Node correlation
		for (int i = 0; i < n.size(); i++)
		{
			nodes[i].n = n[i];
			int node_number = n[i]->number;						// **should** be i

			for (int j = 0; j < e.size(); j++)
			{
				for (int k = 0; k < nnod; k++)
				{
					if ( e[j]->node[k]->number == node_number )
					{
						nodes[i].e.push_back( e[j] );
					}
				}
			}
		}

		cout << "element correlation" << endl;

		// 1. Element correlation
		for (int i = 0; i < e.size(); i++)
		{
			elements[i].e = e[i];
			elements[i].f.resize( nface);
			elements[i].edge.resize( 12);

			// test each boundary for elements, 
			// test the face flag in each element for boundaries.
			// should be able to utilize node data structure here to test.
			// THIS IS 2D SPECIFIC: using 2 nodes, node 0 == face. Need to genericize. 
			if (nnod == 4)	// 2D
			{
				face_association_2d(i);
			}
			else if ( nnod == 8)	// 3D
			{
				face_association_3d(i);
				edge_association_3d(i);
			}

				// add faces from element
			for (int j = 0; j < e[i]->face.size(); j++)		// should be ok 2d or 3d if we provide same fn, will sort by input
			{
				if (nnod == 4)	// 2D
				{
					int face = get_face(e[i]->face[j]->n);
					elements[i].f[face].boundary = true;
				}
				else if ( nnod == 8)	// 3D
				{
					int face = get_face_3d(e[i]->face[j]->n);
					elements[i].f[face].boundary = true;
				}
			}
		}
	}

		// appears to work
	void face_association_3d(int i)
	{
		for (int j = 0; j < nface; j++)		// face iter
		{
			Tensor<int, 1> n = get_nodes_from_face_3d(j);

			// if we have element match(es) in nodes j and k, add them to the face.
			// may be 0 elements (a boundary face), 1 (same refine level) or more than 1 (higher refined ele)

			int node0 = elements[i].e->node[n(0)]->number;
			int node1 = elements[i].e->node[n(1)]->number;
			int node2 = elements[i].e->node[n(2)]->number;
			int node3 = elements[i].e->node[n(3)]->number;

			// iterate over both nodes, find common elements that are not element i.
			for (int k = 0; k < nodes[node0].e.size(); k++)								// node0
			{
				for (int l = 0; l < nodes[node1].e.size(); l++)							// node1
				{
					for (int m = 0; m < nodes[node2].e.size(); m++)						// node2
					{
						for (int n = 0; n < nodes[node3].e.size(); n++)					// node3
						{
							if (   (nodes[node0].e[k]->number == nodes[node1].e[l]->number)		// 1 == 2 == 3 == 4 && != i
								&& (nodes[node1].e[l]->number == nodes[node2].e[m]->number) 
								&& (nodes[node2].e[m]->number == nodes[node3].e[n]->number) 
								&& (nodes[node1].e[l]->number != i) )
							{
							//	cout << "adding element " << nodes[node0].e[k]->number << " to face " << j << " of element " << i << endl;
								elements[i].f[j].e.push_back( nodes[node0].e[k]);
							}
						}
					}
				}
			}
		}	
	//	cin.get();
	}

		// same phenomena does not occur in 2d although we may be interested in it.
	void edge_association_3d(int i)
	{
		// create array of edge data
		std::vector< std::vector< int > > edge;
		edge.resize( 12);
		edge[0] += 0, 1;		// bottom edge
		edge[1] += 1, 2;
		edge[2] += 2, 3;
		edge[3] += 3, 0;

		edge[4] += 0, 4;		// center edge
		edge[5] += 1, 5;
		edge[6] += 2, 6;
		edge[7] += 3, 7;

		edge[8] += 4, 5;		// top edge
		edge[9] += 5, 6;
		edge[10] += 6, 7;
		edge[11] += 7, 4;
		
		for (int j = 0; j < edge.size(); j++)		// face iter
		{
			Tensor<int, 1> n(2);
			n(0) = edge[j][0];
			n(1) = edge[j][1];

			// if we have element match(es) in nodes j and k, add them to the face.
			// may be 0 elements (a boundary face), 1 (same refine level) or more than 1 (higher refined ele)

			int node0 = elements[i].e->node[n(0)]->number;
			int node1 = elements[i].e->node[n(1)]->number;

			// iterate over both nodes, find common elements that are not element i.
			// and not face members
			for (int k = 0; k < nodes[node0].e.size(); k++)								// node0
			{
				for (int l = 0; l < nodes[node1].e.size(); l++)							// node1
				{
					if (   (nodes[node0].e[k]->number == nodes[node1].e[l]->number)		// 1 == 2 == 3 == 4 && != i
						&& (nodes[node1].e[l]->number != i) )
					{
						int m = nodes[node0].e[k]->number;
						bool match = false;
						for (int n = 0; n < elements[i].f.size(); n++)
						{
							for (int ne = 0; ne < elements[i].f[n].e.size(); ne++)
							{
								if (elements[i].f[n].e[ne]->number == m)
									match = true;
						//		cout << "face: " << elements[i].f[n].e[ne]->number << " match? " << match << endl;
							}
						}
						if (!match)
						{
						//	cout << "adding edge " << m << " for element " << i << endl;
							elements[i].edge[j].e.push_back( nodes[node0].e[k]);
						}
					}
				}
			}
		}	
//		cin.get();
	
	}

	void face_association_2d(int i)
	{
		for (int j = 0; j < nface; j++)		// face iter
		{
			Tensor<int, 1> n = get_nodes_from_face(j);

			// if we have element match(es) in nodes j and k, add them to the face.
			// may be 0 elements (a boundary face), 1 (same refine level) or more than 1 (higher refined ele)

			int node0 = elements[i].e->node[n(0)]->number;
			int node1 = elements[i].e->node[n(1)]->number;

			// iterate over both nodes, find common elements that are not element i.
			for (int k = 0; k < nodes[node0].e.size(); k++)
			{
				for (int l = 0; l < nodes[node1].e.size(); l++)
				{
					if ( (nodes[node0].e[k]->number == nodes[node1].e[l]->number) && (nodes[node1].e[l]->number != i) )
					{
							elements[i].f[j].e.push_back( nodes[node0].e[k]);
					}
				}
			}
		}	
	}

	int nnod;
	int nface;

	std::vector< ele<Elements> > elements;
	std::vector< node<Elements, Nodes> > nodes;
};


#pragma once

#include <iostream>
using namespace std;

#include "../FDV/Function/CalcLength.h"
#include "ElementAssociation.h"

	// Finite Element
#include "../FiniteElement/get_face.h"
#include "../FiniteElement/TestFunctions.h"				// test

#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

struct PriorRefinePair
{
	int elno;					// element bordering the refine element
	int face;					// face on the refine element
	std::vector< int> node;		// nodes on the bordering element -> try converting to Tensor, figure out what the problems are.
};								// don't do a fullscale transition try making a temp struct, assigning, check and compare ...

/*
		if (iter == 0 || iter == 2 || iter == 4 || iter == 6)		// 4 iterations works just fine.
		{
				// Clear matrices in each element
			for_each(elements.begin(), elements.end(), ClearElement<Element *>());

				// Attempt mesh refine
			std::vector<int> elelist;			// elements to be refined
			std::vector<int> refine_level;		// breakpoints for each level refine
			AdaptiveRefineList<Element,Node>(elelist, refine_level, elements, nodes, nnod, nface);
			MeshRefine2D<Element, Node>(elements, nodes, elelist, refine_level, ndim, neqn, nnod,iter, nface);
			cout << "Done adapting mesh" << endl;
		}
*/
/*		if (iter == 10)		/// mesh unrefine is broken
		{
				// Attempt mesh unrefine
			MeshUnrefine2D<Element, Node>(elements, nodes, ndim, neqn, nnod);

			cout << "done unrefining mesh" << endl;
			cin.get();
		}
*/

	/// Refine a 2D/3D Mesh
template<class Elements, class Nodes> struct MeshRefine2D
{
	MeshRefine2D(std::vector<Elements*> & elements, std::vector<Nodes*> & nodes, 
		std::vector<int> & elelist, std::vector<int> & refine_level, int ndim, int neqn, int nnod, int iter, int nface)
		: elements(elements),
		  nodes(nodes),
		  ndim(ndim),
		  neqn(neqn),
		  nnod(nnod),
		  nface(nface),
  		  ea(elements, nodes, nnod, nface)
	{
		if (nnod == 4)
		{
			nrnod = 5;		// number of refined nodes
		}
		else if (nnod == 8)
		{
			nrnod = 19;
		}

		int n = refine_level.size() - 2;
		cout << "n: " << n << endl;

		for (int j = 0; j <= n; j++)
		{
			// might need to re-map element relations between 0 and 1 evals ... 
			for (int i = refine_level[j]; i < refine_level[j+1]; i++)			// works so long as we append elements, not insert.
			{
				cout << "refining " << i << " of " << elelist.size() << " (" << elements[ elelist[i] ]->number << ")" << endl;
				refine(elements[ elelist[i] ], elelist, n);
			}

			if (j != n)
				ea.reset(elements, nodes);		// update element associations with each refine phase.
		}
	}

	void refine(Elements * e, std::vector<int> & refinelist, int max_refine_level)
	{
		// 0. Determine what nodes to add by determining if neighbors are being refined...
		// nodelist: Boolean vector of faces bordering this element that match criteria (need refining)
		Tensor<bool, 1> nodelist = findNeighbors(elements, e->number, refinelist, max_refine_level);

		cout << "nodelist: " << nodelist << endl;

			// If there are two refined elements, elelist will return -1.
			// Go in and do a deeper search checking for 2-element faces.
		std::vector<PriorRefinePair> priorRefine = addPriorRefine(nodelist, e);

		cout << "w/refine: " << nodelist << endl;

			/// Add in boundary refinements, if any
		addBoundaries(nodelist, e);

		cout << "w/bc:     " << nodelist << endl;

		// 1. Add new nodes, add to global node stack
		//    assign node numbers to array "new_node" for elements to use when making Hnm / node decisions
		// If we had a vector of zeta components we could make this a single function that loops over zetas,
		// using nodelist as a truth variable. In 3D we'd have to pre-parse nodelist. Something to think about
		// to reduce visible code size and repetition.

		Tensor<double, 1> new_node_number(nrnod);
		if (nnod == 4)
		{
			create_nodes_2d(nodelist, e, new_node_number);
		}
		else if (nnod == 8)
		{
			create_nodes_3d(nodelist, e, new_node_number);
		}

		cout << "final:    " << nodelist << endl;

			// create elements
		createElements(e, nodelist, new_node_number);

			// eventually this could be condensed into one function if we use a struct
			// could intelligently handle the stuff. For now function it away.
		if (nnod == 4)
		{
			do_hanging_nodes_2d(e, new_node_number, priorRefine);
		}
		else if (nnod == 8)
		{
			do_hanging_nodes_3d(e, new_node_number, priorRefine);
		}
	};

	void do_hanging_nodes_3d(Elements * e, Tensor<double, 1> & new_node, std::vector<PriorRefinePair> & priorRefine
		)
	{
		for (int i=0; i < priorRefine.size(); i++)
		{
			cout << "updating element " << priorRefine[i].elno << " with face data " << endl;

			std::vector< int> local = priorRefine[i].node;
			int face = priorRefine[i].face;
			int elno = priorRefine[i].elno;

			/// NEW METHOD WORKS - residuals match to 5 iters, 2 refines.
			/// 1. Figure out if we share a center node using all 4 nodes. (0.25)
			for (int N = 0; N < 8; N++)
			{
				if (
					    elements[elno]->Hnm(N, local[0]) == 0.25 
					 && elements[elno]->Hnm(N, local[1]) == 0.25
					 && elements[elno]->Hnm(N, local[2]) == 0.25 
					 && elements[elno]->Hnm(N, local[3]) == 0.25
					)
				{
					if ( nodes[new_node(face)]->number == 0 ||  nodes[new_node(face)]->number == -1 )
					{
						cout << "Prescribing 0 or -1 as a new node!" << endl;
						cout << "face: " << face << endl;
					}
					else
					{
						update_hanging_nodes(elements[elno], nodes[new_node(face)], N);		// verify face.		
					}
				}
			}

			/// 2. Iterate over all 4 nodes in pairs, and determine if we share an edge (0.50)
			for (int N = 0; N < 8; N++)
			{
				for (int ni=0; ni<4; ni++)
				{
					for (int nj=ni+1; nj<4; nj++)
					{
						if (
								elements[elno]->Hnm(N, local[ni]) == 0.50 
							 && elements[elno]->Hnm(N, local[nj]) == 0.50
							)
						{
							// Can't just reflect nodes because node face ordering is right hand rule.
							// but most faces (except top / bottom) show exploitable symmetry.
							// can probably fix, but until then, here's a hack: if 0,2 increment, if 1,3 decrement
							// Adjusted face 5 in get_face. need to pick a standard and go from there.
							// may be able to find a better way than this, but it'll probably be a lookup table.

							Tensor<int, 1> n = get_nodes_from_face_3d( get_opposing_face_3d( get_face_3d(local) ) );

							int ni0, nj0;
							if ( ni == 0 || ni == 2)
							{
								ni0 = ni+1;
							}
							else
							{							
								ni0 = ni-1;
							}
							if ( nj == 0 || nj == 2)
							{
								nj0 = nj+1;
							}
							else
							{							
								nj0 = nj-1;
							}

							int edge = get_edge_from_nodes( n(ni0), n(nj0) );
							
							if ( nodes[new_node(edge)]->number == 0 ||  nodes[new_node(edge)]->number == -1 )
							{
								cout << "Prescribing 0 or -1 as a new node!" << endl;
								cout << " edge:            " << edge << endl;
								cout << " using nodes:     " << ni0 << " " << nj0 << endl;
								cout << " translated from: " << ni << " " << nj << endl;
								cout << " local:           " << local[0] << " " << local[1] << " " << local[2] << " " << local[3] << endl;
								cout << " opposing face:   " << n << endl;
							}
							else
							{
								update_hanging_nodes(elements[elno], nodes[new_node(edge)], N);				
							}
						}				
					}
				}
			}
		}	
	}

	void do_hanging_nodes_2d(Elements * e, Tensor<double, 1> & new_node, std::vector<PriorRefinePair> & priorRefine	)
	{
		for (int i=0; i < priorRefine.size(); i++)
		{
			int node_number = 0;

			std::vector< int> local = priorRefine[i].node;
			int elno = priorRefine[i].elno;
			int face = priorRefine[i].face;

			node_number = new_node( face );

			if (elements[elno]->refine_level == e->refine_level)	// same because im now refined too
			{	
				int indexi = -1;
				for (int j=0; j < nnod; j++)
				{
					if ( elements[elno]->Hnm(j,local[0]) == 0.5 && elements[elno]->Hnm(j,local[1]) == 0.5 )
						indexi = j;
				}
				if ( indexi == -1)
				{
					cout << "did not find a match for indexi in Hnm in element " << e->number << endl;
					cout << elements[elno]->Hnm << endl;
					cin.get();
					return;
				}

				update_hanging_nodes(elements[elno], nodes[node_number], indexi);
			}
		}	
	}

		/// Update a hanging node in an adapted element.
		// Fix Hnm, if [I] set adap = false.
		// Adds new node to node[] in element.
		// may need to do more than one face at some point ... like leading edge of a flat plate ... think about it.
	void update_hanging_nodes(Elements * e, Nodes * node, int indexi)
	{
		for (int i=0; i < nnod; i++)
		{
			e->Hnm(indexi, i) = 0;
		}
		e->Hnm(indexi, indexi) = 1.0;

			// Update adap flag if Identity matrix
		Tensor<double, 2> ident = Identity<double>(nnod);      // gcc doesn't like instantiation in if()
		if ( e->Hnm == ident )
		{
			e->adap = false;
		}

			// update element node pointers
		e->node[indexi] = node;
	}

	void create_nodes_2d(Tensor<bool, 1> & nodelist, Elements * e, Tensor<double, 1> & new_node)
	{
		int size = nodes.size();

		Tensor<double,1> zeta(2);
		if ( nodelist(0) ) { zeta(0) =  0.0; zeta(1) = -1.0; new_node(0) = createNodeRef(zeta, e); }
		if ( nodelist(1) ) { zeta(0) =  1.0; zeta(1) =  0.0; new_node(1) = createNodeRef(zeta, e); }
		if ( nodelist(2) ) { zeta(0) =  0.0; zeta(1) =  1.0; new_node(2) = createNodeRef(zeta, e); }
		if ( nodelist(3) ) { zeta(0) = -1.0; zeta(1) =  0.0; new_node(3) = createNodeRef(zeta, e); }

			// and then, regardless, create the central node
		nodelist(4) = true;
		zeta(0) = 0; zeta(1) = 0; new_node(4) = createNodeRef(zeta, e);
	}

		/// Create the nodes in 3D space per the nodelist, and update nodelist for face edge nodes
	void create_nodes_3d(Tensor<bool, 1> & nodelist, Elements * e, Tensor<double, 1> & new_node)
	{
		int size = nodes.size();

		Tensor<double,1> zeta(3);

			// Faces (6)
		if ( nodelist(0) ) { zeta(0) =  0.0; zeta(1) =  0.0; zeta(2) = -1.0; new_node(0) = createNodeRef(zeta, e); }
		if ( nodelist(1) ) { zeta(0) =  0.0; zeta(1) = -1.0; zeta(2) =  0.0; new_node(1) = createNodeRef(zeta, e); }
		if ( nodelist(2) ) { zeta(0) =  1.0; zeta(1) =  0.0; zeta(2) =  0.0; new_node(2) = createNodeRef(zeta, e); }
		if ( nodelist(3) ) { zeta(0) =  0.0; zeta(1) =  1.0; zeta(2) =  0.0; new_node(3) = createNodeRef(zeta, e); }
		if ( nodelist(4) ) { zeta(0) = -1.0; zeta(1) =  0.0; zeta(2) =  0.0; new_node(4) = createNodeRef(zeta, e); }
		if ( nodelist(5) ) { zeta(0) =  0.0; zeta(1) =  0.0; zeta(2) =  1.0; new_node(5) = createNodeRef(zeta, e); }

			// Interior Node (1)
		nodelist(6) = true;
		zeta(0) = 0; zeta(1) = 0; zeta(2) = 0; new_node(6) = createNodeRef(zeta, e);
	
			// Vertices (12)
			// along face 0
		if ( nodelist(0) && nodelist(1) ) { nodelist( 7) = true; zeta(0) =  0.0; zeta(1) = -1.0; zeta(2) = -1.0; new_node(7)  = createNodeRef(zeta, e); }
		if ( nodelist(0) && nodelist(2) ) { nodelist( 8) = true; zeta(0) =  1.0; zeta(1) =  0.0; zeta(2) = -1.0; new_node(8)  = createNodeRef(zeta, e); }
		if ( nodelist(0) && nodelist(3) ) { nodelist( 9) = true; zeta(0) =  0.0; zeta(1) =  1.0; zeta(2) = -1.0; new_node(9)  = createNodeRef(zeta, e); }
		if ( nodelist(0) && nodelist(4) ) { nodelist(10) = true; zeta(0) = -1.0; zeta(1) =  0.0; zeta(2) = -1.0; new_node(10) = createNodeRef(zeta, e); }

			// between middle faces
		if ( nodelist(4) && nodelist(1) ) { nodelist(11) = true; zeta(0) = -1.0; zeta(1) = -1.0; zeta(2) =  0.0; new_node(11) = createNodeRef(zeta, e); }
		if ( nodelist(1) && nodelist(2) ) { nodelist(12) = true; zeta(0) =  1.0; zeta(1) = -1.0; zeta(2) =  0.0; new_node(12) = createNodeRef(zeta, e); }
		if ( nodelist(2) && nodelist(3) ) { nodelist(13) = true; zeta(0) =  1.0; zeta(1) =  1.0; zeta(2) =  0.0; new_node(13) = createNodeRef(zeta, e); }
		if ( nodelist(3) && nodelist(4) ) { nodelist(14) = true; zeta(0) = -1.0; zeta(1) =  1.0; zeta(2) =  0.0; new_node(14) = createNodeRef(zeta, e); }

			// along face 5
		if ( nodelist(5) && nodelist(1) ) { nodelist(15) = true; zeta(0) =  0.0; zeta(1) = -1.0; zeta(2) =  1.0; new_node(15) = createNodeRef(zeta, e); }
		if ( nodelist(5) && nodelist(2) ) { nodelist(16) = true; zeta(0) =  1.0; zeta(1) =  0.0; zeta(2) =  1.0; new_node(16) = createNodeRef(zeta, e); }
		if ( nodelist(5) && nodelist(3) ) { nodelist(17) = true; zeta(0) =  0.0; zeta(1) =  1.0; zeta(2) =  1.0; new_node(17) = createNodeRef(zeta, e); }
		if ( nodelist(5) && nodelist(4) ) { nodelist(18) = true; zeta(0) = -1.0; zeta(1) =  0.0; zeta(2) =  1.0; new_node(18) = createNodeRef(zeta, e); }

	}

		///	New version preemptively checks for existing node
		///	If so, return the existing node
		/// If not, push back new node and return pointer
	int createNodeRef( Tensor<double, 1> & zeta, Elements * e)
	{
		Nodes * node = createNode(zeta, e);

		//double tol = 10E-6;										/// tolerence for matching repeat nodes
		double tol = 10E-4;										/// tolerence for matching repeat nodes

		std::vector<int> discard_node;
		for (int i=0; i < nodes.size(); i++)					// node 1 iterator
		{
//			if(	rms(nodes[i]->x, node->x) < tol )				// FIXME rms works on gcc but not vs. Make tensor impl
			double distance = 0;
			distance = sqrt((nodes[i]->x(0) - node->x(0))*(nodes[i]->x(0) - node->x(0))
				+ (nodes[i]->x(1) - node->x(1))*(nodes[i]->x(1) - node->x(1)));
			if(distance < tol)
			{
				delete node;									// UNVERIFIED
				return nodes[i]->number;
			}
		}
			/// If we get here, no match, return the node we made.
		node->number = nodes.size();
		nodes.push_back( node);
		return node->number;
	}

		/// Create a new node by interpolation of the original element
	Nodes * createNode( Tensor<double, 1> & zeta, Elements * e)
	{
		Nodes * node = new Nodes(ndim, neqn);

		Tensor<double, 1> Phi(nnod);
		if ( nnod == 4)
		{
			Phi = testfunction(zeta);
		}
		if ( nnod == 8)
		{
			Phi = testfunction3d(zeta);
		}

		node->x.clear();
		for (int N = 0; N < nnod; N++)
		{
			for (int i=0; i < neqn; i++)
			{
				node->U(i) += Phi(N) * e->node[N]->U(i);
				node->U0(i) += Phi(N) * e->node[N]->U0(i);
			}
			for (int i=0; i < ndim; i++)
			{
				node->x(i) += Phi(N) * e->node[N]->x(i);			
			}
		}

		// Determine Dirichlet data
		addDirichlet(zeta, e, node);

		return node;
	}

		/// Create elements, call specific 2D / 3D function for nodal correlations.
		/// Do the generic work here.
		/// booleans for set_element_node and Hnm... figure it out later.
	void createElements(Element * e, Tensor<bool, 1> & nodelist, Tensor<double, 1> & new_node)
	{
		std::vector<Element *> element;		// need to delete new members

		for (int i=0; i < nnod-1; i++)
		{
			element.push_back(new Elements(nnod, neqn, ndim));
			element[i]->number = elements.size() + i;
			element[i]->Hnm = Identity<double>(nnod);
		}
		e->Hnm = Identity<double>(nnod);

		// set all nodes in all ele's =e
		for (int i = 0; i < element.size(); i++)
		{
			for (int j = 0; j < nnod; j++)
			{
				element[i]->node[j] = e->node[j];
			}
		}

			/// Do specific stuff
		if (nnod == 4)
		{
			createElements_2d(e, element, nodelist, new_node);
		}
		else if (nnod == 8)
		{
			createElements_3d(e, element, nodelist, new_node);
		}

			// Calculate new area/length in element
		CalcLength<Elements *> calc_length(ndim,4);				// nnod hardcoded but that's ok we only use quads
		for (int i=0; i < element.size(); i++)
		{
			calc_length(element[i]);
		}
		calc_length(e);

			// flag elements as refined
		e->refine_level++;
		for (int i=0; i < nnod-1; i++)
		{
			element[i]->refine_level = e->refine_level;
		}

			// push back onto stack
		for (int i=0; i < element.size(); i++)
		{
			elements.push_back(element[i]);
		}
	}

		/// Create three new elements given the knowledge of the system
		/// can probably make a generic createElements function? all is the same except for the 
		/// booleans for set_element_node and Hnm... figure it out later.
	void createElements_2d(Element * e, std::vector<Element *> & element, Tensor<bool, 1> & nodelist, Tensor<double, 1> & new_node)
	{
		// create nodelist - Hnm analogies
		std::vector< std::vector<double> > Hnm;
		Hnm.resize( 5);

		// using analogies derived Friday 10/24 on paper, create Hnm rows that correspond to nodelist numbers
		////////	  0	    1     2     3 
		Hnm[ 0] += 0.50, 0.50, 0.00, 0.00; 		// 1. interp face 0
		Hnm[ 1] += 0.00, 0.50, 0.50, 0.00; 		// 2. interp face 1
		Hnm[ 2] += 0.00, 0.00, 0.50, 0.50; 		// 3. interp face 2
		Hnm[ 3] += 0.50, 0.00, 0.00, 0.50; 		// 4. interp face 3
		Hnm[ 4] += 0.00, 0.00, 0.00, 0.00; 		// center node SHOULD NOT BE USED

				/// ELEMENT 0
		//					element   local node, new node, blargh
		set_Element_Node(element[0], 0,  0, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 2,  1, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 3,  4, Hnm, nodelist, new_node);

			/// ELEMENT 1
		set_Element_Node(element[1], 1,  1, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 3,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 0,  4, Hnm, nodelist, new_node);

			/// ELEMENT 2
		set_Element_Node(element[2], 0,  3, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 2,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 1,  4, Hnm, nodelist, new_node);

			/// ELEMENT e
		set_Element_Node(e, 1,  0, Hnm, nodelist, new_node);
		set_Element_Node(e, 3,  3, Hnm, nodelist, new_node);
		set_Element_Node(e, 2,  4, Hnm, nodelist, new_node);

			// Add face data to the element
		addFaces_2d(e, element);
	}

		/// Create seven new elements given the knowledge of the system (3D)
	void createElements_3d(Element * e, std::vector<Element *> & element, Tensor<bool, 1> & nodelist, Tensor<double, 1> & new_node)
	{
		// create nodelist - Hnm analogies
		std::vector< std::vector<double> > Hnm;
		Hnm.resize( 19);

		// using analogies derived Friday 10/24 on paper, create Hnm rows that correspond to nodelist numbers
		////////	 0	   1     2     3     4     5     6     7
		Hnm[ 0] += 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00; 		// 1. interp face 0
		Hnm[ 1] += 0.25, 0.25, 0.00, 0.00, 0.25, 0.25, 0.00, 0.00; 		// 1. interp face 1
		Hnm[ 2] += 0.00, 0.25, 0.25, 0.00, 0.00, 0.25, 0.25, 0.00; 		// 1. interp face 2
		Hnm[ 3] += 0.00, 0.00, 0.25, 0.25, 0.00, 0.00, 0.25, 0.25; 		// 1. interp face 3
		Hnm[ 4] += 0.25, 0.00, 0.00, 0.25, 0.25, 0.00, 0.00, 0.25; 		// 1. interp face 4
		Hnm[ 5] += 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25; 		// 1. interp face 5
		Hnm[ 6] += 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00; 		// 1. Middle node - NEVER INTERPOLATED
		Hnm[ 7] += 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00; 		// 1. Vertice between face 0,1
		Hnm[ 8] += 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00; 		// 1. Vertice between face 0,2
		Hnm[ 9] += 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00; 		// 1. Vertice between face 0,3
		Hnm[10] += 0.50, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00; 		// 1. Vertice between face 0,4
		Hnm[11] += 0.50, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00; 		// 1. Vertice between face 4,1
		Hnm[12] += 0.00, 0.50, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00; 		// 1. Vertice between face 1,2
		Hnm[13] += 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.50, 0.00; 		// 1. Vertice between face 2,3
		Hnm[14] += 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.50; 		// 1. Vertice between face 3,7
		Hnm[15] += 0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00; 		// 1. Vertice between face 5,1
		Hnm[16] += 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.00; 		// 1. Vertice between face 5,2
		Hnm[17] += 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.50; 		// 1. Vertice between face 5,3
		Hnm[18] += 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.50; 		// 1. Vertice between face 5,4

		/// ELEMENT 0
		//					element   local node, new node, blargh
		set_Element_Node(element[0], 0,  7, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 2,  8, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 3,  0, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 4,  1, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 5, 12, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 6,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[0], 7,  6, Hnm, nodelist, new_node);

		/// ELEMENT 1
		set_Element_Node(element[1], 0,  0, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 1,  8, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 3,  9, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 5,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 6, 13, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 7,  3, Hnm, nodelist, new_node);
		set_Element_Node(element[1], 4,  6, Hnm, nodelist, new_node);

		/// ELEMENT 2
		set_Element_Node(element[2], 0, 10, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 1,  0, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 2,  9, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 4,  4, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 6,  3, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 7, 14, Hnm, nodelist, new_node);
		set_Element_Node(element[2], 5,  6, Hnm, nodelist, new_node);

		/// ELEMENT 3
		set_Element_Node(element[3], 0, 11, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 1,  1, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 3,  4, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 5, 15, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 6,  5, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 7, 18, Hnm, nodelist, new_node);
		set_Element_Node(element[3], 2,  6, Hnm, nodelist, new_node);

		/// ELEMENT 4
		set_Element_Node(element[4], 0,  1, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 1, 12, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 2,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 4, 15, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 6, 16, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 7,  5, Hnm, nodelist, new_node);
		set_Element_Node(element[4], 3,  6, Hnm, nodelist, new_node);

		/// ELEMENT 5
		set_Element_Node(element[5], 1,  2, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 2, 13, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 3,  3, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 4,  5, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 5, 16, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 7, 17, Hnm, nodelist, new_node);
		set_Element_Node(element[5], 0,  6, Hnm, nodelist, new_node);

		/// ELEMENT 6
		set_Element_Node(element[6], 0,  4, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 2,  3, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 3, 14, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 4, 18, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 5,  5, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 6, 17, Hnm, nodelist, new_node);
		set_Element_Node(element[6], 1,  6, Hnm, nodelist, new_node);

		/// ELEMENT e
		set_Element_Node(e, 1,  7, Hnm, nodelist, new_node);
		set_Element_Node(e, 2,  0, Hnm, nodelist, new_node);
		set_Element_Node(e, 3, 10, Hnm, nodelist, new_node);
		set_Element_Node(e, 4, 11, Hnm, nodelist, new_node);
		set_Element_Node(e, 5,  1, Hnm, nodelist, new_node);
		set_Element_Node(e, 7,  4, Hnm, nodelist, new_node);
		set_Element_Node(e, 6,  6, Hnm, nodelist, new_node);

			// Add face data to the element
		addFaces_3d(e, element);
	}

		/// Generic function to set the node in an element
		/// If it is a member of the nodelist, set the new node.
		/// If not, keep the old (original element) node and interpolate using provided Hnm
	void set_Element_Node(Element * e, int elenode, int node, std::vector< std::vector<double> > & Hnm, Tensor<bool, 1> & nodelist, Tensor<double, 1> & new_node)
	{
		if( nodelist(node) )
		{
			e->node[elenode] = nodes[ new_node(node) ];			
		}
		else
		{
				// Hnm matrix manipulation
			e->adap = true;
			for (int i = 0; i < nnod; i++)
			{
				e->Hnm(elenode, i) = Hnm[node][i];
			}
		}
	}

		/// Calculate dirichlet boundary conditions via quadrature.
		/// Nodes lying directly between 2 dirichlet nodes evaluate to one, all others are zero'd.
		/// @param zeta the isoparametric interpolation coordinates
		/// @param e the element being interpolated
		/// @param n the node where the results will be calculated
	void addDirichlet( Tensor<double, 1> & zeta, Elements * e, Nodes * n)
	{
		Tensor<double, 1> dirichlet(neqn);
		Tensor<double, 1> Phi(nnod);
		if ( nnod == 4)
		{
			Phi = testfunction(zeta); 
		}
		else if (nnod == 8)
		{
			Phi = testfunction3d(zeta);
		}

		/// Need to sum into a double, otherwise the partial values all evaluate (true) == 1.
		for (int N = 0; N < nnod; N++)
		{
			for (int i=0; i < neqn; i++)
			{
				dirichlet(i) += Phi(N) * (double)e->node[N]->dirichlet(i);
			}
		}

			// Convert from double to boolean value.
		for (int i=0; i < neqn; i++)
		{
			if (dirichlet(i) < 1)
			{
				n->dirichlet(i) = false;
			}
			else
			{
				n->dirichlet(i) = true;
			}
		}
	}

		/// This is 2D spefic
		/// Add faces to new elements, extrapolating from the original element
		/// Remove unaplicable faces from the original element.
		/// NOTE: this is wrong, deleting inline like this bumps i+1 down and you skip it ...
		/// Change to be like 3D which I think is correct ... 
	void addFaces_2d(Element * e, std::vector<Element *> & element)
	{
		if (e->number == 29)
		{
			cout << "modifying element 29" << endl;
		}
//		for (int i = 0; i < e->face.size(); i++)
		for (int i = e->face.size()-1; i>=0; i--)
		{
			int face = get_face(e->face[i]->n);
			if ( face == 0 )		// Bottom face, implement in elements e and 0
			{ 
				element[0]->face.push_back( e->face[i] );
			}
			
			else if ( face == 1 )	// Right face, implement in elements 0 and 1 and remove from e
			{
				element[0]->face.push_back( e->face[i] );
				element[1]->face.push_back( e->face[i] );
				e->face.erase(e->face.begin()+i, e->face.begin()+i+1);
			}
			
			else if ( face == 2 )	// Top face, implement in elements 1 and 2 and remove from e
			{
				element[1]->face.push_back( e->face[i] );
				element[2]->face.push_back( e->face[i] );
				e->face.erase(e->face.begin()+i, e->face.begin()+i+1);
			}
			
			else if ( face == 3 )	// Left face, implement in elements 3 and e
			{
				element[2]->face.push_back( e->face[i] );
			}			
		}
	}

	/// Should be OK need to validate by running something of substance to completion ... 
	/// fixed a glaring problem (deleted faces inline, means you skipped faces)
	void addFaces_3d(Element * e, std::vector<Element *> & element)
	{
		std::vector<Face *> faces = e->face;
		e->face.clear();

		for (int i = faces.size()-1; i >= 0; i--)
		{
			int face = get_face_3d(faces[i]->n);
			if ( face == 0 )							// Bottom face, implement in elements e and 0,1,2
			{ 
				e->face.push_back( faces[i] );
				element[0]->face.push_back( faces[i] );
				element[1]->face.push_back( faces[i] );
				element[2]->face.push_back( faces[i] );
			}
			else if ( face == 4 /*1*/ )						// Front face, implement in elements e and 2,3,6
			{
				e->face.push_back( faces[i] );
				element[2]->face.push_back( faces[i] );
				element[3]->face.push_back( faces[i] );
				element[6]->face.push_back( faces[i] );
			}
			else if ( face == 1 /*2*/ )						// Right face, implement in elements e and 1,4,5
			{
				e->face.push_back( faces[i] );
				element[0]->face.push_back( faces[i] );
				element[3]->face.push_back( faces[i] );
				element[4]->face.push_back( faces[i] );
			}
			else if ( face == 2 /*3*/ )						// Back face, implement in elements 0,1,4,5 and remove from e
			{
				element[0]->face.push_back( faces[i] );
				element[1]->face.push_back( faces[i] );
				element[4]->face.push_back( faces[i] );
				element[5]->face.push_back( faces[i] );
			}			
			else if ( face == 3 /*4*/ )						// Left face, implement in elements 1,2,5,6 and remove from e
			{
				element[1]->face.push_back( faces[i] );
				element[2]->face.push_back( faces[i] );
				element[5]->face.push_back( faces[i] );
				element[6]->face.push_back( faces[i] );
			}
			else if ( face == 5 )						// Top face, implement in elements 3,4,5,6 and remove from e
			{
				element[3]->face.push_back( faces[i] );
				element[4]->face.push_back( faces[i] );
				element[5]->face.push_back( faces[i] );
				element[6]->face.push_back( faces[i] );
			}
			else
			{
				cout << "face not matched: " << i << endl;
			}
		}
	}

		/// Inspect the element for faces, if any, to determine if the element lies on a boundary. 
		/// If so, add those nodes to the nodelist as they can be refined.
		/// @param nodelist list of nodal positions to be refines
		/// @param e element being refined
	void addBoundaries(Tensor<bool, 1> & nodelist, Elements * e)
	{
		for (int i = 0; i < e->face.size(); i++)
		{
			if (nnod == 4)
			{
				nodelist( get_face( e->face[i]->n ) ) = true;
			}
			if (nnod == 8)
			{
				nodelist( get_face_3d( e->face[i]->n ) ) = true;
			}
		}
	}

		/// Determine elements bordering faces that will be at the same refine_level once we are refined.
		/// These faces will then share nodes with us, so we need to update our correlation (nodelist)
	std::vector<PriorRefinePair> addPriorRefine(Tensor<bool, 1> & nodelist, Element * e)
	{
		std::vector<PriorRefinePair> priorRefine;

		int i = e->number;

		int bordering_eles = 0;
		if (nnod == 4) bordering_eles = 2;
		else if (nnod == 8) bordering_eles = 4;

		for (int j = 0; j < ea.elements[i].f.size(); j++)
		{
			if (ea.elements[i].f[j].e.size() == bordering_eles)		// 2 for 2D, 4 for 3D
			{
				nodelist(j) = true;

				for (int k = 0; k < ea.elements[i].f[j].e.size(); k++)
				{
					Tensor<int, 1> nodes(nnod / 2);

					if ( nnod == 4)
						nodes = get_nodes_from_face( get_opposing_face(j) );		// 2D
					else if ( nnod == 8)
						nodes = get_nodes_from_face_3d( get_opposing_face_3d(j) );	// 3D
	
					vector<int> nodesvec;
					for (int l = 0; l < nodes.imax(); l++)
					{
						nodesvec += nodes(l);
					}
					
					PriorRefinePair pair = {ea.elements[i].f[j].e[k]->number, j, nodesvec};

					priorRefine.push_back(pair);
				}
			}
		}

		return priorRefine;
	}

		///	Finds and returns nearby elements that remain to be refined.
		/// @param e element vector
		/// @param i element of interest in the vector
		/// @param facelist retun Tensor of elements by true/false
	Tensor<bool, 1> findNeighbors(std::vector<Elements *> & e , int i, std::vector<int> & refinelist, int max_refine_level)
	{
		Tensor<bool, 1> face(nrnod);
		face.clear();

		int jlist = 0;
		for (int j=0; j < refinelist.size(); j++)
		{
			if (refinelist[j] == i)
			{
				jlist = j;
			}
		}
		for (int j = 0; j < ea.elements[i].f.size(); j++)
		{
			for (int ki = 0; ki < ea.elements[i].f[j].e.size(); ki++)
			{
				int k = ea.elements[i].f[j].e[ki]->number;

				for (int l = 0; l < refinelist.size(); l++)			// should be able to simplify this for/if/if/else expression
				{
					if (refinelist[l] == k)
					{
						if ((e[k]->refine_level == e[i]->refine_level) && (l < jlist))
						{}
						else
						{
							face(j) = true;
						}
					}
				}
			}
		}

		return face;
	}

	int nrnod;									/// number of refined nodes, ie, the number of nodes (potentially) created by refinement (5 or 19, 2D, 3D respectively)
	int ndim, neqn, nnod, nface;

	std::vector<Nodes*> & nodes;
	std::vector<Elements*> & elements;
	ElementAssociation<Element, Node> ea;		/// Element association, correlates adjacent elements and nodes
};

// namespace Adaptive
// {
		/// Handles modification of test functions in an adaptive mesh.
	void adap(Tensor<double, 1> & Phi, Tensor<double, 2> & dPxi, Tensor<double, 2> & Hnm)
	{
		Tensor<double, 1> PhiHnm(Phi.imax());
		Tensor<double, 2> dPxiHnm(dPxi.imax(), dPxi.jmax());

		PhiHnm.clear();
		dPxiHnm.clear();

		for (int N = 0; N < Phi.imax(); N++)
		{
			for (int M = 0; M < Phi.imax(); M++)
			{
				PhiHnm(M)  += Phi(N) * Hnm(N, M);

				for (int j = 0; j < dPxi.jmax(); j++)
				{
					dPxiHnm(M, j) += dPxi(N, j) * Hnm(N, M);
				}
			}
		}
	
		Phi  = PhiHnm;
		dPxi = dPxiHnm;
	}

		// adapt a 3D face (2D surface) from 3D element (3D nodes)
	void adap_face(Tensor<double, 1> & phibou, Tensor<double, 2> & dpidx, std::vector<Node *> & node, Element * e)
	{
		// phi = phi_n * H_nm
		// but only over the nodes of interest.
		// need a int tensor flag.
		Tensor<int, 2> correlation(4,2);
		correlation.clear();

		// load correlation
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (node[j]->number == e->node[i]->number)
				{
					for (int k = 0; k < 3; k++)
					{
						correlation(j, 0) = j;
						correlation(j, 1) = i;
					}
				}
			}
		}

		// now, multiply, using correlation as a truth test. 
		Tensor<double, 1> phi_t( phibou.imax());	
		Tensor<double, 2> dphi_t( dpidx.imax() , dpidx.jmax() );
		phi_t.clear();
		dphi_t.clear();

		for (int N = 0; N < 4; N++)
		{
			for (int M = 0; M < 4; M++)
			{
				int Nlocal  = correlation(N, 0);
				int Mlocal  = correlation(M, 0);
				int Nglobal = correlation(N, 1);
				int Mglobal = correlation(M, 1);

				phi_t(Mlocal)  += phibou(Nlocal) * e->Hnm(Nglobal, Mglobal);
				for (int j = 0; j < 2; j++)		// BUG, SHOULD BE =3
				{
					dphi_t(Mlocal, j) += dpidx(Nlocal, j) * e->Hnm(Nglobal, Mglobal);
				}
			}
		}
			
		phibou  = phi_t;
		dpidx = dphi_t;		
	}
// };


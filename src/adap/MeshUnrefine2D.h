#pragma once

#include <iostream>
#include <algorithm>
#include "CalcLength.h"
using namespace std;

//#include <boost/assign/std/vector.hpp>	// for 'operator+=()'
#include "boost/assign/std/vector.hpp"
using namespace boost;
using namespace boost::assign;			// bring 'operator+=()' into scope

static bool sort_using_greater_than(int u, int v)
{
   return u > v;
}

static bool sort_using_less_than(int u, int v)
{
   return u < v;
}

/*		if (iter == 10)		/// mesh unrefine is broken
{
// Attempt mesh unrefine
MeshUnrefine2D<Element, Node>(elements, nodes, ndim, neqn, nnod);

cout << "done unrefining mesh" << endl;
cin.get();
}
*/
	/// Unrefine a 2D Mesh
template<class Elements, class Nodes> 
struct MeshUnrefine2D
{
	ElementAssociation<Element, Node> e;		/// Element association, correlates adjacent elements and nodes

	MeshUnrefine2D(std::vector<Elements*> & elements, std::vector<Nodes*> & nodes, int ndim, int neqn, int nnod)
		: elements(elements),
		  nodes(nodes),
		  ndim(ndim),
		  neqn(neqn),
		  nnod(nnod),
		  e(elements, nodes, nnod, 4)
	{
		// determine highest refine_level
		int size = elements.size();
		int n = 0;
		for (int i=0; i < size; i++)
		{
			if (elements[i]->refine_level > n)
				n = elements[i]->refine_level;
		}

		// find elements with a refine level > 1 and a s1 < 0.2 (0.6 for testing)
		std::vector<int> refine_list;
		for (int i=0; i < size; i++)
		{
			if (elements[i]->refine_level > 1 && elements[i]->s1 < 0.6)
			{
				refine_list.push_back( i);
			}
		}

		// find if we have any quadrants
		// should do it by node, I think ... 
		Tensor<int, 1> ele(4);
		for (int i = 0; i < 4; i++)
			ele(i) = -1;
		for (int i = 0; i < nodes.size(); i++)
		{
			if (e.nodes[i].e.size() == nnod)				// don't refine along a boundary, won't work and algorithm breaks.
			{
				for (int j = 0; j < nnod; j++)
				{
					for (int k = 0; k < refine_list.size(); k++)
					{
						if (e.nodes[i].e[j]->number == e.elements[k].e->number )
						{
							ele(j) = k;
						}
					}
				}
				bool match = true;
				for (int j = 0; j < ele.imax; j++)
				{
					if ( ele(j) == -1)
						match = false;
				}
				if(match)
				{
					cout << "Match found " << ele << endl;
					unrefine( ele );
				}
			}
		}

//		discard_nodes();

		discard_elements();
	}

		// This logic works but watch out if orderings or things change in MeshRefine.
	void unrefine(Tensor<int, 1> & ele)
	{
		element_unrefine(ele);
	}

	void element_unrefine(Tensor<int, 1> & ele)
	{
		// we want to keep the first element and discard the last three <--- DO NOT NECESSARILY KNOW THIS unless we sort low->high
		// and to be safe should geometrically verify.

		int i = ele(0);

		// 1. reset node1, node2, node3 in element 0.
		elements[i]->node[1] = elements[ ele(1) ]->node[1];
		elements[i]->node[2] = elements[ ele(2) ]->node[2];
		elements[i]->node[3] = elements[ ele(3) ]->node[3];

		// 2. re-calculate length of first element
		CalcLength<Elements *> calc_length(ndim);
		calc_length(elements[i]);

		// 3. aggregate faces
		for (int j = 1; j < ele.imax; j++)
		{
			for (int k = 0; k < elements[ ele(j) ]->face.size(); k++)
			{
				int bc = elements[ ele(j) ]->face[k].bc;
				bool bcflag = false;
				for (int l = 0; l < elements[i]->face.size(); l++)
				{
					if (elements[i]->face[l].bc == bc)
						bcflag = true;
				}
				if (!bcflag)
				{
					elements[i]->face.push_back( elements[ ele(j) ]->face[k] );
				}
			}
		}

		// 4. clear Hnm, adap, refine_level
		elements[i]->Hnm.clear();
		elements[i]->adap = false;
		elements[i]->refine_level--;

		// 5. mark elements for deletion
		for (int j = 1; j < ele.imax; j++)
		{
			discard_element.push_back( ele(j) );
		}

		// 6. mark nodes for deletion - need to be careful because Hnm matrices could
		//    be pointing to a corner node here. Deferred to later.
		std::vector< int> nodes; 
		nodes += elements[ ele(0) ]->node[1]->number,
				 elements[ ele(0) ]->node[2]->number,
				 elements[ ele(0) ]->node[3]->number,
			     elements[ ele(1) ]->node[2]->number,
			     elements[ ele(2) ]->node[3]->number;


		 for (int j = 0; j < nnod; j++)
		 {
			for (int k = 0; k < nodes.size(); k++)
			{
				if ( elements[i]->node[j]->number == nodes[k] )
				{
					cout << "node " << nodes[k] << " should not be deleted" << endl;
					nodes[k] = -1;
				}
			}
		 }
		for (int j = 0; j < nodes.size(); j++)
		{
			if ( nodes[j] != -1)
			{
		 		discard_node.push_back( nodes[j] );
			}
		}

	}

		// All the nodes deleted here are past the original domain
		// So we don't run into any problems, etc
	void discard_nodes()
	{
		std::sort(discard_node.begin(), discard_node.end(), sort_using_greater_than);
		
		for (int i=0; i < discard_node.size(); i++)
		{
			cout << "deleting node " << discard_node[i] << endl;
			nodes.erase(nodes.begin()+discard_node[i], nodes.begin()+discard_node[i]+1);
		}
	}

	void discard_elements()
	{
		std::sort(discard_element.begin(), discard_element.end(), sort_using_greater_than);

		// since they are sorted from highest to lowest, removing an element should not invalidate the list.

		for (int i=0; i < discard_element.size(); i++)
		{
			elements.erase(elements.begin()+discard_element[i], elements.begin()+discard_element[i]+1);
		}
	}

	int ndim, neqn, nnod;

	std::vector<int> discard_node;
	std::vector<int> discard_element;

	std::vector<Nodes*> & nodes;
	std::vector<Elements*> & elements;
};


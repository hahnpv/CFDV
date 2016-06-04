#pragma once

// init CFD stuff with MPI.

// presently 2D specific, to be extended to 3D.

#include <vector>
#include "Utility/dictionary.h"
#include <iostream>
using namespace std;
struct MPI_init_CFD
{
		int nodemin;
		int nodemax;
		int neqn;
	MPI_init_CFD(int argc, char * argv[], std::vector<Element *> & elements, std::vector<Node *> & nodes, int neqn, int nnod, dictionary & key)
		: elements(elements),	// for functions
		  nnod(nnod),
		  neqn(neqn),
		  key(key),
		  nodes(nodes)
	{
		//// MPI Initialize ////
		int rc = MPI_Init(&argc,&argv);
		size = MPI::COMM_WORLD.Get_size();
		rank = MPI::COMM_WORLD.Get_rank();
/*
		int min = 0;
		int max = 0;
		get_elements(min, max);

		nodemin = 0;
		nodemax = 0;
		get_nodes(nodemin, nodemax);

		get_gmres(neqn, nodemin);

		stat();
		*/
}

	void stat()
	{
		cout << "MPI Init CFD    " << endl;
		cout << "Rank            " << rank << " of " << size << endl;
		cout << "Total Elements: " << elements.size() << endl;
		cout << "From:           " << elements[0]->number << " to " << elements[elements.size()-1]->number << endl;
		cout << "Total Nodes:    " << nodesize / neqn << endl;
		cout << "From:           " << nodemin << " to " << nodemax << endl;
		cout << "Offset:         " << offset << endl;
//		cout << "Node Iterator:  " << nodes[0]->number << " to " << nodes[/*0*/ nodesize/neqn]->number /*+ nodesize / neqn*/ << endl;
		cout << endl;	
	}

	void get_elements(int & min, int & max)
	{
			// Chop elements
		int window = key.get_val<int>("elements") / size;
		int elemin = window * rank;

		// Instead of chopping, we'd then load these elements.
		if (rank != size-1)		// prevent roundoff errors from excluding elements in last rank
		{
			elements.erase( elements.begin() + elemin + window, elements.end() );	// right
		}

		elements.erase( elements.begin(), elements.begin() + elemin );				// left

		// binary load (elemin, window); (start, size)
	}

	void get_nodes(/*int & nodemin, int & nodemax*/)
	{
			// This search is fine, if we aggregate it from the elements we are loading.
		// Keep nodes relevant to elements
		nodemin = elements[0]->node[0]->number;
		nodemax = elements[0]->node[0]->number;
		for (int i = 0; i < elements.size(); i++)
		{
			for (int j = 0; j < nnod; j++)
			{
				if ( elements[i]->node[j]->number < nodemin)
					nodemin = elements[i]->node[j]->number;
				else if ( elements[i]->node[j]->number > nodemax)
					nodemax = elements[i]->node[j]->number;
			}
		}
		nodes.erase(nodes.begin() + nodemax+1, nodes.end());
		nodes.erase(nodes.begin(), nodes.begin() + nodemin);

	//	cout << "Chopping to nodes " << nodemin << " " << nodemax << " in rank " << rank << endl;

		// binary load (nodemin, nodemax)
	}

	void get_gmres(/*int neqn, int nodemin*/)
	{
			//// GMRES Partitioning ////
		nodesize  = 0;

		nedgenode_left = 0;
		nedgenode = 0;

		int rank_left_max  = 0;
		int rank_right_min = 0;

		//// SEND BLOCK ////
		if (rank != 0)
		{
			cout << "rank " << rank << " sending to rank " << rank-1 << endl;
			MPI::COMM_WORLD.Send( &nodes[0]->number, 1, MPI_INT, rank-1, 1);					// rank_right_min
		}
		if (rank != size-1)
		{
			cout << "rank " << rank << " receiving from rank " << rank+1 << endl;
			MPI::COMM_WORLD.Send( &nodes[nodes.size()-1]->number, 1, MPI_INT, rank+1,  1);		// rank_left_max
		}
		//// SEND BLOCK ////

		//// RECV BLOCK ////
		if (rank != size-1)
		{
			MPI::COMM_WORLD.Recv( &rank_right_min, 1, MPI_INT, rank+1, 1);
		}
		if (rank != 0)
			MPI::COMM_WORLD.Recv( &rank_left_max, 1, MPI_INT, rank-1, 1);
		//// RECV BLOCK ////

		int ne_l = 0;
		int ne_r = 0;
		if ( rank != 0)
		{
			ne_l = rank_left_max - nodes[0]->number;
		}
		if ( rank != size-1)
		{
			ne_r = nodes[nodes.size()-1]->number - rank_right_min;
		}
		cout << "diff left: " << ne_l << " right: " << ne_r << endl;

		if (rank != 0)						// add in the missing node ... 
			nedgenode_left = ne_l + 1;
		if (rank != size-1)
			nedgenode      = ne_r + 1;

		cout << "nedgenode (left, right): " << nedgenode_left << " " << nedgenode << endl;

		// this works, 2D
		if (rank < size-1)
		{
			nodesize = (nodes.size() - nedgenode) * neqn;
		}
		else
		{
			nodesize = nodes.size() * neqn;
		}

		offset = nodemin;	// duh, take the number of nodes chopped ... 

		// find the nodes to iterate over for RMSError.
		NodeIteratorStart = nodes.begin();

		if (rank == size-1)
			NodeIteratorEnd = nodes.end();// - nedgenode;
		else		// if rank == 0
			NodeIteratorEnd = nodes.end() - nedgenode;
	}

private:

	dictionary & key;
	std::vector<Node *> & nodes;
	std::vector<Element *> & elements;

public:
	int nodesize;				/// nodes * eqns, try to replace with just nodes and work out in programs.
	int rank;					/// our rank number
	int size;					/// total number of ranks
	int offset;					/// offset in node count (GMRES)
	int nedgenode;				/// number of nodes on the edge bordering the next rank
	int nedgenode_left;
	int nnod;
	std::vector<Node *>::iterator NodeIteratorStart;
	std::vector<Node *>::iterator NodeIteratorEnd;
};
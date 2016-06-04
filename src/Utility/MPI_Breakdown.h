#pragma once

// init CFD stuff with MPI.

// 2D specific

// Break down the element / node / GMRES points so that it doesn't need to be calc'd on the fly
// will allow for more efficient loading, etc.


#include <vector>
// #include "dictionary.h"
#include <iostream>
using namespace std;
struct MPI_Breakdown
{
	static void clear(std::string path)
	{
		ofstream fout((path + "//MPIBreakdown").c_str(),ios::out);
		fout.close();
	}
	std::vector<int> nodemin;
	std::vector<int> nodemax;
	std::vector<int> min;
	std::vector<int> max;
	std::string path;
	int neqn;
	MPI_Breakdown( std::vector<Element *> & elements, std::vector<Node *> & nodes, int neqn, int nnod, std::string path, int size)
		: elements(elements),	// for functions
		  nnod(nnod),
		  neqn(neqn),
		  key(key),
		  nodes(nodes),
		  path(path),
		  size(size)
	{
		nodemin.resize( size);
		nodemax.resize( size);
		min.resize( size);
		max.resize( size);

		cout << "calculating elements and nodes" << endl;
		for (int rank = 0; rank < size; rank++)
		{
			get_elements(min[rank], max[rank], rank);
			get_nodes(min[rank], max[rank], nodemin[rank], nodemax[rank]);
		}
		std::vector<int> ne_l;
		std::vector<int> ne_r;
		std::vector<int> nodesize;
		std::vector<int> offset;
		std::vector<int> NodeIteratorStart;
		std::vector<int> NodeIteratorEnd;

		ne_l.resize( size);
		ne_r.resize( size);
		nodesize.resize( size);
		offset.resize( size);
		NodeIteratorStart.resize( size);
		NodeIteratorEnd.resize( size);

		cout << "Calculating GMRES" << endl;
		get_gmres(ne_l, ne_r, nodesize, offset,	NodeIteratorStart, NodeIteratorEnd);

		cout << "stats: " << endl;
		for (int rank = 0; rank < size; rank++)
		{
			cout << "Rank " << rank << endl;
			cout << "Number of elements: " << max[rank] - min[rank] << endl;
			cout << "From:               " << min[rank] << " " << max[rank] << endl;
			cout << "Nodes from:         " << nodemin[rank] << " " << nodemax[rank] << endl;
			cout << "Nedgenode:          " << ne_l[rank] << " " << ne_r[rank] << endl;
			cout << "NodeIterator:       " << nodemin[rank] << " " << nodemax[rank] - ne_r[rank] << endl;
		}

		ofstream fout( (path + "//MPIBreakdown").c_str(),ios::app);
		fout << size << endl;
		fout << "elemin, elemax, nodemin, nodemax, nedgenode_l, nedgenode_r, NodeIteratorMin, NodeIteratorMax" << endl;
		for (int rank = 0; rank < size; rank++)
		{
			fout << min[rank] << " " << max[rank] << " " << nodemin[rank] << " " << nodemax[rank] << " " << ne_l[rank] << " " << ne_r[rank] << " " << nodemin[rank] << " " << nodemax[rank] - ne_r[rank] << endl;
		}

		fout.close();
	}

	void get_elements(int & min, int & max, int rank)
	{
			// Chop elements
		int window = /*key.get_val<int>("elements")*/  elements.size() / size;
		int elemin = window * rank;

		min = elemin;
		max = elemin + window;		// max is one past size (ie: delete criteria)
		if ( rank == size-1)
			max = elements.size();
	}

	void get_nodes(int elemin, int elemax, int & min, int & max)
	{
			// This search is fine, if we aggregate it from the elements we are loading.
		// Keep nodes relevant to elements
		min = elements[elemin]->node[0]->number;
		max = elements[elemin]->node[0]->number;
		for (int i = elemin; i < elemax; i++)
		{
			for (int j = 0; j < nnod; j++)
			{
				if ( elements[i]->node[j]->number < min)
					min = elements[i]->node[j]->number;
				else if ( elements[i]->node[j]->number > max)
					max = elements[i]->node[j]->number;
			}
		}
	}

	void get_gmres(std::vector<int> & ne_l, std::vector<int> & ne_r, std::vector<int> & nodesize, std::vector<int> & offset,
		std::vector<int> & NodeIteratorStart, std::vector<int> & NodeIteratorEnd)
	{
			//// GMRES Partitioning ////
		for (int i = 0; i < size; i++)
		{
			if ( i != size-1)
				ne_r[i] = nodemax[i] - nodemin[i+1] + 1;
			if ( i != 0)
				ne_l[i] = nodemax[i-1] - nodemin[i] + 1;
		}

		for (int i = 0; i < size; i++)
		{
			if ( i != size-1)
			{
				nodesize[i] = (nodemax[i] - nodemin[i] - ne_r[i]) * neqn;
			}
			else
			{
				nodesize[i] = (nodemax[i] - nodemin[i]) * neqn;
			}

			offset[i] = nodemin[i];
		}

		// find the nodes to iterate over for RMSError.
		for (int i = 0; i < size; i++)
		{
			NodeIteratorStart[i] = nodemin[i];
			NodeIteratorEnd[i]   = nodemax[i] - ne_r[i];
		}
	}

private:

	dictionary & key;
	std::vector<Node *> & nodes;
	std::vector<Element *> & elements;

public:
	int nodesize;				/// nodes * eqns, try to replace with just nodes and work out in programs.
//	int rank;					/// our rank number
	int size;					/// total number of ranks
	int offset;					/// offset in node count (GMRES)
//	int nedgenode;				/// number of nodes on the edge bordering the next rank
//	int nedgenode_left;
	int nnod;
//	std::vector<Node *>::iterator NodeIteratorStart;
//	std::vector<Node *>::iterator NodeIteratorEnd;
};
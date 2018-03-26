#pragma once

// init CFD stuff with MPI.

// 2D specific

// Break down the element / node / GMRES points so that it doesn't need to be calc'd on the fly
// will allow for more efficient loading, etc.


#include <vector>
#include <iostream>
using namespace std;
struct MPI_Breakdown
{
	static void clear(std::string path)
	{
		ofstream fout((path + "//MPIBreakdown").c_str(),ios::out);
		fout.close();
	}
	
	// TEST
	static void read(std::string path, int rank, int & size, int & nodemin, int & nodemax, int & elemin, int & elemax, int & nedgenode, int & nedgenode_left, int & max, int & node_iterator_min, int & node_iterator_max)
	{
		// still not as parallel as it could be but much better.
		cout << "path: " << path << endl;
		ifstream mbd((path + "//MPIBreakdown").c_str(), ios::in);
		std::vector<int> data;
		for (int i = 1; i <= 16; i++)		// temp hardcode, doesnt work except for 2?
		{
			std::string line;
			getline(mbd, line);		// config #
			cout << rank << " config: " << line << endl;
			getline(mbd, line);		// description
			cout << rank << " description: " << line << endl;
			for (int j = 0; j < i; j++)	// read in breakdown points
			{
				getline(mbd, line);
				cout << rank << " dat: " << line << endl;
				data = split<int>(line, " ");

				if (i == size && j == rank)
				{
					cout << "rank = " << rank << " using line: " << line << endl;
					break;

					cout << rank << " data: ";
					for (unsigned int i = 0; i < data.size(); i++)
						cout << data[i] << " ";
					cout << endl;
				}
			}
			if (i == size)
				break;
		}
		mbd.close();

		elemin = data[0];
		elemax = data[1];

		nodemin = data[2];
		nodemax = data[3];

		nedgenode = data[5];
		nedgenode_left = data[4];
		max = (data[7] - data[6]) + 1;

		node_iterator_min = data[6];
		node_iterator_max = data[7];
	}
	
	std::vector<unsigned int> nodemin;
	std::vector<unsigned int> nodemax;
	std::vector<int> min;
	std::vector<int> max;
	std::string path;
	int neqn;
	MPI_Breakdown( std::vector<Element *> & elements, std::vector<Node *> & nodes, int neqn, int nnod, std::string path, int size)
		: elements(elements),	// for functions
		  nnod(nnod),
		  neqn(neqn),
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
		int window = elements.size() / size;
		int elemin = window * rank;

		min = elemin;
		max = elemin + window;		// max is one past size (ie: delete criteria)
		if ( rank == size-1)
			max = elements.size();
	}

	void get_nodes(int elemin, int elemax, unsigned int & min, unsigned int & max)
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
	std::vector<Node *> & nodes;
	std::vector<Element *> & elements;

public:
	int nodesize;				/// nodes * eqns, try to replace with just nodes and work out in programs.
	int size;					/// total number of ranks
	int offset;					/// offset in node count (GMRES)
	int nnod;
};
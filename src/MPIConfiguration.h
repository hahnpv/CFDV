#pragma once

#include "ConfigurationBase.h"
#include "MPI_Init_CFD.h"
#include "Utility/IO/MPI_TecplotOut.h"
#include "FDV/Function/MPI_RMSError.h"
#include "Utility/Solvers/MPI_GMRES.h"
//#include "MPI_Save.h"

struct MPIConfiguration : ConfigurationBase
{
	MPIConfiguration(int argc, char * argv[], std::vector<Element *> & elements, std::vector<Node *> & nodes, dictionary & key)
		: ConfigurationBase(key)
	{
		int    iter    = key.get_val<int>("iter");
		double cfl     = key.get_val<double>("cfl");					/// Prescribed CFL number
		int itermax    = key.get_val<int>("itermax");					/// Maximum number of time iterations
		int gmres_iter = key.get_val<int>("gmresiter");				/// Number of GMRES iterations 
		int gmres_rest = key.get_val<int>("gmresrestart");				/// Number of GMRES restarts

		MPI_init_CFD initMPI(argc, argv, elements, nodes, neqn, nnod, key);
//MPI_init_CFD initMPI(argc, argv, elements, nodes, neqn, nnod, key);

		// tweak initMPI (or move parts here) to chop element/node, but instead, use it to load element/node.
		// LoadBinaryData init(args, elements, nodes);							/// new binary method, change constructor to not autoconstruct
		int min = 0;
		int max = 0;
		initMPI.get_elements(min, max);

		initMPI.get_nodes();

		initMPI.get_gmres();

		initMPI.stat();

//		init.read_adap();

		rank = initMPI.rank;
		size = initMPI.size;
		offset = initMPI.offset;
		nedgenode = initMPI.nedgenode;
		nedgenode_left = initMPI.nedgenode_left;
		int nodesize = initMPI.nodesize;

		NodeIteratorStart = initMPI.NodeIteratorStart;
		NodeIteratorEnd = initMPI.NodeIteratorEnd;

		gmres = new MPI_GMRES( gmres_rest, gmres_iter, nodesize, neqn, nnod, ndim, size, rank, offset);
	}
	
	void Output(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, double time, bool ele_data)
	{
		MPI_TecplotOut(elements, nodes, iter, ndim, time, false, rank,offset,size, /*nedgenode,*/ nedgenode_left, nnod, neqn);
	}

	void RMSErr(double t, int iter)
	{
		MPI_RMSError<Node *> rmserr("RMSError.txt", t, iter, neqn);
		for_each(NodeIteratorStart, NodeIteratorEnd, std::tr1::ref(rmserr));
	}

	void Save(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, dictionary & key)
	{
	//	MPI_Save(elements, nodes, ndim, nnod, iter, rank, key);
	}
	
	int size;
	int offset;
	int nedgenode;
	int nedgenode_left;
};
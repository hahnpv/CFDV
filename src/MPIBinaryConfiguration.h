#pragma once

#include "ConfigurationBase.h"
#include "MPI_Init_CFD.h"
#include "Utility/IO/MPI_TecplotOut.h"
#include "FDV/Function/MPI_RMSError.h"
#include "Utility/Solvers/MPI_GMRES.h"
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
#include "Utility/IO/ASCIItoBinary/SaveBinary.h"
//#include "MPI_Save.h"

#include "Utility/MPI_Breakdown.h"

#include <iomanip>
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#include "boost/filesystem/fstream.hpp"   // filestreams, makes avail to other classes and may cause problems?

#include <functional> // std::ref

struct MPIBinaryConfiguration : ConfigurationBase
{
	MPIBinaryConfiguration(int argc, char * argv[], std::vector<Element *> & elements, std::vector<Node *> & nodes, dictionary & key, LoadBinaryData & init)
		: ConfigurationBase(key)
	{
		int    iter    = key.get_val<int>("iter");
		double cfl     = key.get_val<double>("cfl");					/// Prescribed CFL number
		int itermax    = key.get_val<int>("itermax");					/// Maximum number of time iterations
		gmres_iter = key.get_val<int>("gmresiter");					/// Number of GMRES iterations 
		gmres_rest = key.get_val<int>("gmresrestart");				/// Number of GMRES restarts

		int rc = MPI_Init(&argc,&argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		// instead of MPI_init_CFD, read in the MPIBreakdown file...
		// still not as parallel as it could be but much better.
		path = key.get_val<std::string>("path");
		cout << "path: " << path << endl;
		ifstream mbd( (path + "//MPIBreakdown").c_str(), ios::in);
		std::vector<int> data;
		for (int i = 1; i <= 8; i++)		// temp hardcode, doesnt work except for 2?
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
			if ( i == size)
				break;
		}
		mbd.close();

		elemin = data[0];
		elemax = data[1];
		
		nodemin = data[2];
		nodemax = data[3];

		init.read_nodes(nodemin, nodemax);		
		init.read_elements(elemin, elemax);	

		cout << "nodes.size: " << nodes.size() << endl;

		nedgenode      = data[5];
		nedgenode_left = data[4];
		NodeIteratorStart = nodes.begin();
		NodeIteratorEnd = nodes.begin() + (data[7] - data[6]) + 1;
		max = (data[7] - data[6]) + 1;

		offset = nodemin;
		nodesize = (nodes.size() - nedgenode) * neqn;
		init.read_adap();

		cout << "rank:          " << rank << endl;
		cout << "element range: " << elemin << " " << elemax << " " << elements.size() << endl;
		cout << "node range:    " << nodemin << " " << nodemax << " " << nodes.size() << endl;
		cout << "nedgenode:     " << nedgenode << " left: " << nedgenode_left << endl;
		cout << "node iter:     " << 0 << " to " << data[7] - data[6] + 1 << endl;//" which is " 
		cout << "offset:        " << offset << endl;
		cout << "node size:     " << nodesize << endl;
							//	 << nodes[0]->number << " to " << nodes[(data[7] - data[6]) + 1]->number << endl;

		gmres = new MPI_GMRES( gmres_rest, gmres_iter, nodesize, neqn, nnod, ndim, size, rank, offset);
	}
	

	void reset_gmres()
	{
		delete gmres;
		gmres = new MPI_GMRES( gmres_rest, gmres_iter, nodesize, neqn, nnod, ndim, size, rank, offset);
	}

	void Output(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, double time, bool ele_data)
	{
		MPI_TecplotOut(elements, nodes, iter, ndim, time, ele_data , rank, offset, size, /*nedgenode,*/ nedgenode_left, nnod, neqn);
	}

	void RMSErr(double t, int iter)
	{
		MPI_RMSError<Node *> rmserr("RMSError.txt", t, iter, neqn);
		for_each(NodeIteratorStart, NodeIteratorEnd, ref(rmserr));
	}

	void Save(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, dictionary & key)
	{
		std::string savepath = path + "//" + to_string<int>(iter) + "//";
		boost::filesystem::create_directory(savepath);
		std::string fname = savepath + "binarynodes";
		// need to save using MPI functionale
		// just need to save nodes - element, face correlations remains the same, EXCEPT when you grid refine!
		// eventually dump an updated control file too

		// TODO FIXME note that these node_t, ele_t are defined in LoadBinary, use those!

		///////////////////
		// Nodes -> 2d node
		///////////////////
		MPI_Datatype TType;		/// node data xfer

		if ( ndim == 3)
			MPI_Type_contiguous(8, MPI_DOUBLE, &TType);
		else if (ndim == 2)
			MPI_Type_contiguous(6, MPI_DOUBLE, &TType);
		else if (ndim == 1)
			MPI_Type_contiguous(4, MPI_DOUBLE, &TType);
		MPI_Type_commit(&TType);

		struct node_t			// TODO check i think we need to make this x[ndim]
		{
			double x[2];
			double rho;
			double v[2];
			double T;
		};

		cout << "writing nodes to Save file " << rank << endl;
		int nele = nodes.size(); // 10 + rank;					// nodes.size(), nele should be nnod but this var is taken
		cout << "writing  " << max << " nodes "  << nele << endl;
//		nele = max;
		node_t * buf;
		buf = new node_t[nele];				// elements are unique

		for (int i = 0; i < nele; i++)		// get data from nodes
		{
			buf[i].rho = nodes[i]->rho;
			buf[i].T   = nodes[i]->T;
			for (int j = 0; j < ndim; j++)
			{
				buf[i].x[j] = nodes[i]->x(j);
				buf[i].v[j] = nodes[i]->v(j);
			}
		}

		cout << "opening file" << endl;
		int size = 0;
		MPI_Allreduce(&nodemax, &size, 1,MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,
			fname.c_str(),  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_File_set_size(fh, size * sizeof(node_t) );

		cout << "writing shared" << endl;
		MPI_File_seek(fh, sizeof(node_t)* nodes[0]->number /*nele*rank*/, MPI_SEEK_SET);			// naive as we have overlapping nodes use node iter
		MPI_File_write(fh, buf, nele, TType, MPI_STATUS_IGNORE); 

		cout << "closing" << endl;
		MPI_File_close(&fh);
		cout << "done" << endl;

		/////////////////
		// elements
		/////////////////
		//std::string 
		
		fname = savepath + "binaryelements";

		if (ndim == 3)
			MPI_Type_contiguous(24, MPI_INT, &TType);
		else if (ndim == 2)
			MPI_Type_contiguous(8, MPI_INT, &TType);
		else if (ndim == 1)
			MPI_Type_contiguous(2, MPI_INT, &TType);
		MPI_Type_commit(&TType);

		struct ele_t
		{
			int node[4];
			int bc[4];
		};
		
		cout << "writing elements to Save file " << rank << endl;
		nele = elements.size(); // 10 + rank;					// nodes.size(), nele should be nnod but this var is taken
		cout << "writing  " << max << " nodes " << nele << endl;
		//nele = max;
		ele_t * ebuf;
		ebuf = new ele_t[nele];				// elements are unique

		for (int i = 0; i < nele; i++)		// get data from nodes
		{
			ebuf[i].bc[0] = 0;
			ebuf[i].bc[1] = 0;
			ebuf[i].bc[2] = 0;
			ebuf[i].bc[3] = 0;

			for (int j = 0; j < nnod; j++)		// nodes
			{
				ebuf[i].node[j] = elements[i]->node[j]->number;
			}
			for (unsigned int j = 0; j < elements[i]->face.size(); j++)		// faces fixme
			{
				ebuf[i].bc[elements[i]->face[j]->face] = elements[i]->face[j]->bc;
			}
		}

		cout << "opening file" << endl;
		size = 0;
		MPI_Allreduce(&elemax, &size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

//		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,
			fname.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_File_set_size(fh, size * sizeof(ele_t));

		cout << "writing shared" << endl;
		MPI_File_seek(fh, sizeof(ele_t)* elements[0]->number /*nele*rank*/, MPI_SEEK_SET);			// naive as we have overlapping nodes use node iter
		MPI_File_write(fh, ebuf, nele, TType, MPI_STATUS_IGNORE);

		cout << "closing" << endl;
		MPI_File_close(&fh);
		cout << "done" << endl;

		// note: this (and grid refinement) is NOT yet mpi capable fixme
		boost::filesystem::ofstream adapfile(savepath + "adap");
		for (int i = 0; i < nele; i++)		// get data from nodes
		{
			if (elements[i]->adap) // are adap flags being set???
			{
				// write file
				adapfile << elements[i]->number << " "
					<< elements[i]->adap << " "
					<< elements[i]->refine_level << " "
					<< "double "						// FIXME double probably was when i was thinking of converting to single, 2/8/8 unclear
					<< ndim << " "
					<< nnod << " "
					<< nnod << " ";
				for (int j = 0; j < nnod*nnod; j++)
				{
					adapfile << elements[i]->Hnm(j) << " ";
				}
				adapfile << endl;

			}
		}
		adapfile.close();

		// not sure if this should work for multiple processes
		for(int i=1; i<9;i++)
			MPI_Breakdown(elements, nodes, neqn, nnod, savepath, i);
	}
	
	int size;
	int offset;
	int nedgenode;
	int nedgenode_left;
	int nodesize;
	int nodemin, nodemax;
	int elemin, elemax;
	std::string path;

	int max;				// int equiv of max iterator

	int gmres_iter, gmres_rest;
};

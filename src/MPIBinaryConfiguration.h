#pragma once

#include "ConfigurationBase.h"
#include "MPI_Init_CFD.h"
#include "Utility/IO/MPI_TecplotOut.h"
#include "FDV/Function/MPI_RMSError.h"
#include "Utility/Solvers/MPI_GMRES.h"
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
#include "Utility/IO/ASCIItoBinary/SaveBinary.h"
//#include "MPI_Save.h"

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
					for (int i = 0; i < data.size(); i++)
						cout << data[i] << " ";
					cout << endl;
				}
			}
			if ( i == size)
				break;
		}
		mbd.close();

		int elemin = data[0];
		int elemax = data[1];
		
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
		MPI_TecplotOut(elements, nodes, iter, ndim, time, false, rank,offset,size, /*nedgenode,*/ nedgenode_left, nnod, neqn);
	}

	void RMSErr(double t, int iter)
	{
		MPI_RMSError<Node *> rmserr("RMSError.txt", t, iter, neqn);
		for_each(NodeIteratorStart, NodeIteratorEnd, std::tr1::ref(rmserr));
	}

	void Save(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, dictionary & key)
	{
		std::string fname = "nodes"  + to_string<int>(iter);
		// need to save using MPI functionale
		// just need to save nodes - element, face correlations remains the same
		// eventually dump an updated control file too

		// Nodes -> 3d node
		MPI_Datatype TType;		/// node data xfer

		if ( ndim == 3)
			MPI_Type_contiguous(8, MPI_DOUBLE, &TType);
		else if (ndim == 2)
			MPI_Type_contiguous(6, MPI_DOUBLE, &TType);
		else if (ndim == 1)
			MPI_Type_contiguous(4, MPI_DOUBLE, &TType);
		MPI_Type_commit(&TType);

		struct node_t
		{
			double x[3];
			double rho;
			double v[3];
			double T;
		};

		cout << "writing nodes to Save file " << rank << endl;
		int nele = nodes.size(); // 10 + rank;					// elements.size()
		cout << "writing  " << max << " nodes "  << nele << endl;
		nele = max;
		node_t * buf;
		buf = new node_t[nele];				// elements are unique

		for (int i = 0; i < nele; i++)		// get data from elements
		{
			buf[i].rho = nodes[i]->rho;
			buf[i].T   = nodes[i]->T;
			for (int j = 0; j < ndim; j++)
			{
				buf[i].x[j] = nodes[i]->x(j);
				buf[i].v[j] = nodes[i]->v(j);
			}
		}
//		for (int i = 0; i < nele; i++)		// get data from elements
//			buf[i].rho = 10*rank + i;

		cout << "opening file" << endl;
		int size = 0;
		MPI_Allreduce(&nodemax, &size, 1,MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,
			fname.c_str(),  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_File_set_size(fh, size * sizeof(node_t) );

		cout << "writing shared" << endl;
		MPI_File_seek(fh, sizeof(node_t)* nodes[0]->number /*nele*rank*/, MPI_SEEK_SET);			// naive as we have overlapping nodes use node iter
		MPI_File_write(fh, buf, nele, TType, MPI_STATUS_IGNORE);									// fixme status object unclear 
//		fh.Write_all(buf, nele, TType);

		cout << "closing" << endl;
		MPI_File_close(&fh);
		cout << "done" << endl;
	}
	
	int size;
	int offset;
	int nedgenode;
	int nedgenode_left;
	int nodesize;
	int nodemin, nodemax;
	std::string path;

	int max;				// int equiv of max iterator

	int gmres_iter, gmres_rest;
};
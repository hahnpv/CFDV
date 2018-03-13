#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#include "boost/filesystem/fstream.hpp"   // filestreams, makes avail to other classes and may cause problems?

using namespace std;

#include "LoadBinary.h"
#include "../../../FiniteElement/get_face.h"

struct SaveBinaryData
{
	SaveBinaryData(
			 std::vector< Element * > & elements, 
			 std::vector< Node *> & nodes,
			 int ndim,
			 int neqn,
			 int nnod,
			 int rank,
			 std::string savepath
			 )
			 : ndim(ndim),
			   neqn(neqn),
			   nnod(nnod),
		       rank(rank),
			   elements(elements),
			   nodes(nodes),
			   savepath(savepath)
	{ }

	void write_elements(int elemin, int elemax)
	{
		MPI_Datatype TType;		/// element

		std::string fname = savepath + "binaryelements";

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
		int nele = elements.size(); // 10 + rank;					// nodes.size(), nele should be nnod but this var is taken
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
		int size = 0;
		MPI_Allreduce(&elemax, &size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		//		MPI_File fh;
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,
			fname.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_File_set_size(fh, size * sizeof(ele_t));

		cout << "writing shared" << endl;
		MPI_File_seek(fh, sizeof(ele_t)* elements[0]->number /*nele*rank*/, MPI_SEEK_SET);			// naive as we have overlapping nodes use node iter
		MPI_File_write(fh, ebuf, nele, TType, MPI_STATUS_IGNORE);

		cout << "closing" << endl;
		MPI_File_close(&fh);
		cout << "done" << endl;
	}

	void write_nodes(int nodemin, int nodemax)
	{
		std::string fname = savepath + "binarynodes";

		MPI_Datatype TType;		/// node data xfer

		if (ndim == 3)
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
		node_t * buf;
		buf = new node_t[nele];				// elements are unique

		for (int i = 0; i < nele; i++)		// get data from nodes
		{
			buf[i].rho = nodes[i]->rho;
			buf[i].T = nodes[i]->T;
			for (int j = 0; j < ndim; j++)
			{
				buf[i].x[j] = nodes[i]->x(j);
				buf[i].v[j] = nodes[i]->v(j);
			}
		}

		cout << "opening file" << endl;
		int size = 0;
		MPI_Allreduce(&nodemax, &size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,
			fname.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		MPI_File_set_size(fh, size * sizeof(node_t));

		cout << "writing shared" << endl;
		MPI_File_seek(fh, sizeof(node_t)* nodes[0]->number /*nele*rank*/, MPI_SEEK_SET);			// naive as we have overlapping nodes use node iter
		MPI_File_write(fh, buf, nele, TType, MPI_STATUS_IGNORE);

		cout << "closing" << endl;
		MPI_File_close(&fh);
		cout << "done" << endl;
	}

	void write_adap()
	{
		int nele = elements.size();
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
	}

	int ndim;
	int neqn;
	int nnod;

	int rank;

	std::vector< Element * > elements;
	std::vector< Node *> nodes;
	std::string savepath;
};

#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

#include "../../FDV/Function/FDVParam.h"
#include "../../FDV/Element.h"

// Output a Tecplot FE data file

struct elem
{
	int number;
	int node[8];	// 8
};

struct fdv
{
	double s1;
	double s2;
	double s3;
	double s4;
};

struct pnode
{
	double number;
	double x;
	double y;	// z
	double z;
	double rho;
	double u;
	double v;	// w 
	double w;
	double T;
	double p;
};

struct MPI_TecplotOut
{
	MPI_TecplotOut(std::vector<Element *> &Elements, std::vector<Node *> &Nodes, int iter, int ndim, double time, bool ele_data, int rank, int offset, int size, int nedgenode, int nnod, int neqn)
	{
		/// Define Datatypes
		MPI_Datatype NodeType;
		MPI_Type_contiguous(10, MPI_DOUBLE, &NodeType);
		MPI_Type_commit(&NodeType);

		MPI_Datatype ElementType;
		MPI_Type_contiguous(9, MPI_INT, &ElementType);
		MPI_Type_commit(&ElementType);

		MPI_Datatype FdvType;
		MPI_Type_contiguous(4, MPI_DOUBLE, &FdvType);
		MPI_Type_commit(&FdvType);

		// do it like MPI_GMRES, make a vector of stuff to send and send it, except in this case, no return vec.
		std::vector< elem> elements;
		std::vector< pnode> nodes;
		std::vector< fdv> ffdv;

		int delta = 0;
		if ( rank > 0)
			delta += nedgenode;				// this does NOT WORK in new methodology, need to use other int generated by nedgenode calc in mpi_init.
									// nedgenode_left, this is nedgenode_right.
		// build vec
		// nodes
		for (unsigned int i = delta; i < Nodes.size(); i++)
		{
			pnode n;
			n.number = Nodes[i]->number;
			n.x = Nodes[i]->x(0);
			n.y = Nodes[i]->x(1);
			if (ndim == 3)
			{
				n.z = Nodes[i]->x(2);
			}
			n.rho = Nodes[i]->rho;
			n.u = Nodes[i]->v(0);
			n.v = Nodes[i]->v(1);
			if (ndim == 3)
			{
				n.w = Nodes[i]->v(2);
			}
			n.T = Nodes[i]->T;
			n.p = Nodes[i]->p;

			nodes.push_back( n);
		}

		// elements
		for (unsigned int i = 0; i < Elements.size(); i++)
		{
			elem e;
			e.number = Elements[i]->number;
			for (int x = 0; x < nnod; x++)
			{
				e.node[x] = Elements[i]->node[x]->number;
			}
			elements.push_back( e);

			if (ele_data)
			{
				// fdv
				double s1 = 0;
				double s2 = 0;
				double s3 = 0;
				double s4 = 0;
				FDVParam<Element *> fdvparam;			/// Calculate the FDV parameters
				fdvparam(Elements[i], s1, s2, s3, s4);

				fdv f;
				f.s1 = s1;
				f.s2 = s2;
				f.s3 = s3;
				f.s4 = s4;
				ffdv.push_back(f);
			}
		}


		// send/recv call
		if (rank > 0)
		{
			int tag = 1;
			int other_rank = 0;

			//// SEND BLOCK ////
			int num_msg_send = elements.size();
			MPI_Send( &num_msg_send, 1, MPI_INT, other_rank, tag, MPI_COMM_WORLD);
			for (int j = 0; j < num_msg_send; j++)
			{
				MPI_Send( &elements[j], 1, ElementType, other_rank, tag, MPI_COMM_WORLD);
			}
			num_msg_send = nodes.size();
			MPI_Send( &num_msg_send, 1, MPI_INT, other_rank, tag, MPI_COMM_WORLD);
			for (int j = 0; j < num_msg_send; j++)
			{
				MPI_Send( &nodes[j], 1, NodeType, other_rank, tag, MPI_COMM_WORLD);
			}

			num_msg_send = ffdv.size();
			MPI_Send(&num_msg_send, 1, MPI_INT, other_rank, tag, MPI_COMM_WORLD);
			for (int j = 0; j < num_msg_send; j++)
			{
				MPI_Send(&ffdv[j], 1, FdvType, other_rank, tag, MPI_COMM_WORLD);
			}
			//// SEND BLOCK ////
		}

		if (rank == 0)
		{
			int tag = 1;
			int size = 0;
			MPI_Comm_size(MPI_COMM_WORLD, &size);
			for (int k = 0; k < size - 1; k++)		// size == 2, minus 1 for 0,1, .
			{
				//// RECV BLOCK ////
				int num_msg_recv = 0;
				MPI_Recv( &num_msg_recv, 1, MPI_INT, k+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int j = 0; j < num_msg_recv; j++)
				{
					elem e;
					MPI_Recv( &e, 1, ElementType, k+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					elements.push_back( e);
				}
				MPI_Recv( &num_msg_recv, 1, MPI_INT, k+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int j = 0; j < num_msg_recv; j++)
				{
					pnode n;
					MPI_Recv( &n, 1, NodeType, k+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					nodes.push_back( n);
				}

				MPI_Recv(&num_msg_recv, 1, MPI_INT, k + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int j = 0; j < num_msg_recv; j++)
				{
					fdv f;
					MPI_Recv(&f, 1, FdvType, k + 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					ffdv.push_back(f);
				}
				//// RECV BLOCK ////
			}
		}
		if (rank == 0)
		{
			// write file
			std::ofstream tout(("tecplot_"+to_string<double>(iter)+".plt").c_str(),ios::out);

			if (ndim == 1)
			{
				for(unsigned int n = 0; n < nodes.size(); n++)
				{
					tout << nodes[n].x << " " << nodes[n].y << " " << nodes[n].rho << " " << nodes[n].u << " " << nodes[n].T << " " << nodes[n].p << endl;
				}
			}
			else
			{
				// Write Tecplot header
			tout << "title = \"Output at time = " << time <<" \"" << endl;
			tout << "variables = \"x\", \"y\",";
			if (ndim == 3)
			{
				tout << " \"z\",";
			}
			tout <<	"\"rho\", \"u\", \"v\",";
			if (ndim == 3)
			{
				tout << " \"w\",";
			}
			tout << "\"T\", \"p\",";
			if (ele_data)
			{
				tout << "\"s1\", \"s2\", \"s3\", \"s4\", ";
			}
			tout << endl;
			tout << "zone E=" << elements.size() << ", N="<< nodes.size();
			if (nnod == 3)
			{
				tout << ", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE";
			}
			else
			{
				tout << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL";
			}
			if (ndim == 3)
			{
			//	tout << " ET=BRICK";
				tout << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK";
			}
			tout << endl;
			if (ele_data)
			{
				if ( ndim == 2)
				{
					tout << "VARLOCATION=([8,9,10,11]=CellCentered)" << endl;
				}
				else if (ndim == 3)
				{
					tout << "VARLOCATION=([10,11,12,13]=CellCentered)" << endl;
				}
			}
			tout << endl;
			
				// Set scientific notation w/7 digits of precision
			tout << scientific << endl;
			tout << setprecision(7) << endl;


			// x, y
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].x << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;
			
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].y << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;


			if (ndim == 3)
			{
				for (unsigned int n=0; n<nodes.size(); n++)
				{
					tout << nodes[n].z << " ";

					if (n%100 == 0 && n != 0) tout << endl;
				}
				tout << endl;
			}

			// rho
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].rho << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

					// u
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].u << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

						// v
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].v << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

			if(ndim == 3)
			{
				for (unsigned int n=0; n<nodes.size(); n++)
				{
					tout << nodes[n].w << " ";

					if (n%100 == 0 && n != 0) tout << endl;
				}
				tout << endl;
			}

			// temperature
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].T << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

			// pressure
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n].p << "\t";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

			// s1-s4
			if (ele_data)
			{
				for (unsigned int f = 0; f < ffdv.size(); f++)
				{
					tout << ffdv[f].s1 << "\t";
					if (f % 100 == 0 && f != 0) tout << endl;
				}
				for (unsigned int f = 0; f < ffdv.size(); f++)
				{
					tout << ffdv[f].s2 << "\t";
					if (f % 100 == 0 && f != 0) tout << endl;
				}
				for (unsigned int f = 0; f < ffdv.size(); f++)
				{
					tout << ffdv[f].s3 << "\t";
					if (f % 100 == 0 && f != 0) tout << endl;
				}
				for (unsigned int f = 0; f < ffdv.size(); f++)
				{
					tout << ffdv[f].s4 << "\t";
					if (f % 100 == 0 && f != 0) tout << endl;
				}
			}

			// correlation
			for (unsigned int e=0; e < elements.size(); e++)
			{
				for(int n=0; n < nnod; n++)
				{
					tout << elements[e].node[n]+1 << " ";
				}
				tout << endl;
			}

			}
			tout.close();
		}
	}
};

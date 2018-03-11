#pragma once

#include "mpi.h"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "../../Utility/Tensor.h"
#include "../../FDV/Element.h"
#include "../../Utility/tostring.h"
#include "boost/assign/std/vector.hpp"
using namespace boost::assign;

struct DataStruct
{
	int indice;
	double data[5];		// sending like this chops exec time in HALF
};

struct DataStructAv
{
	int indice;
	int elno;				// specific to Av (Ax)
	int locnode;			// specific to Av (Ax)
	double data[5];			// neqn
};

struct NodeStruct	
{
	double index;
	double node;
	double data[5];
};

	/// the Generalized Minimial RESidual method for solving
	/// the matrix problem Ax = b; this is accomplished using
	/// an element-by-element (EBE) assembly of A and b.
	/// Needs generalization and testing under harsh conditions.
	/// Needs to be more efficient.
	/// Then extension to 3D.
	
#include "GMRESbase.h"

struct MPI_GMRES : public GMRESbase
{
	/// Constructor, intializes dimensions of arrays
	MPI_GMRES(int r, int m, int size, int neqn, int nnod, int ndim, int numthreads, int rank, int offset)
		: size(size),
		  restart(r), 
		  m(m),
		  m_(m),
		  r(size),
		  v(size,m+1),
		  v_tilde(size,m+1),
		  h(m+1,m),
		  g(m+1),
		  y(m),
		  z(size),
		  Ax(size),
		  du(size),
		  b(size),
		  neqn(neqn),
		  nnod(nnod),
		  ndim(ndim),
		  numthreads(numthreads),
		  rank(rank),
		  prefetchnodes(false)
	{
		cout << "MPI_GMRES, rank" << rank << endl;

		delta = offset;
		min = delta;
		rank_down = rank - 1;
		rank_up = rank + 1;
		max = delta + size / neqn - 1;

		cout << "delta (offset): " << delta << endl;
		cout << "min:            " << min << endl;
		cout << "max:            " << max << endl;
		cout << "rank up/down:   " << rank_up << " " << rank_down << endl;
		cout << "size:           " << size << endl;
		cout << "will assign from min to max" << endl;

		/// MPI Datatypes

		DataStruct ds;
		MPI_Aint displ = (size_t)(&ds.data) - (size_t)(&ds);

		int blocklengths[3] = {1,5,1};									// 3 int's, 1 double, 1 eostruct
		MPI_Datatype types[3] = {  MPI_INT, MPI_DOUBLE, MPI_UB };		// datatypes in order
		MPI_Aint lb;
		MPI_Aint intex;
		MPI_Aint displacements[3];								// memory offsets for struct
		MPI_Type_get_extent(MPI_INT, &lb, &intex);						// memory offset from int's to double, dont need double since its last
		displacements[0] = 0;
		displacements[1] = displ;
		displacements[2] = sizeof( ds );
		MPI_Type_create_struct(3, blocklengths, displacements, types, &DataStructType);
		MPI_Type_commit(&DataStructType);

		// Av type
		DataStructAv dsAv;
		displ = (size_t)(&dsAv.data) - (size_t)(&dsAv);

		int blocklengthss[3] = {3,5,1};									// 3 int's, 1 double, 1 eostruct
		MPI_Datatype typess[3] = {  MPI_INT, MPI_DOUBLE, MPI_UB };		// datatypes in order
		MPI_Aint displacementss[3];										// memory offsets for struct
		MPI_Type_get_extent(MPI_INT, &lb, &intex);						// memory offset from int's to double, dont need double since its last
		displacementss[0] = 0;
		displacementss[1] = displ;
		displacementss[2] = sizeof( ds );
		MPI_Type_create_struct(3, blocklengthss, displacementss, typess, &DataStructAvType);
		MPI_Type_commit(&DataStructAvType);

		MPI_Type_contiguous(7, MPI_DOUBLE, &NodeStructType);
		MPI_Type_commit(&NodeStructType);

		do_once = false;		// have yet to prefetch
	}

	/// Solves the system of equations, given the element correlation and node data
	/// @param elements element data structure, provides input data
	void iterate( std::vector< Element *> & elements )
	{
		MPI_Barrier(MPI_COMM_WORLD);
		cout << "MPI_GMRES::Iterate" << endl;

		if (!do_once)
		{
			prefetch(elements);
			do_once = true;
		}

		du.clear();
		
		compose_rhs(elements);				// clear, then compose b

		double tol = MPImag(b) * 1.0E-9;	// calculate tolerence for convergence (to quit gmres restarts early)
		if (rank == 0)
			cout << " tolerence: " << tol << endl;

		for (int gloop = 0; gloop < restart; gloop++)
		{
			clear();						// clear nonpersistant data 

			sum_Ax(elements,Ax,du);			// sum initial residual
			for (unsigned int n = 0; n < r.imax(); n++)
			{
				r(n) = b(n) - Ax(n);
			}

			double mag_r = MPImag(r);		// calculate RMS of residual
			g(0) = mag_r;					// intialize Hessenberg matrix solution vector

			for (int n = 0; n < size; n++)
			{
				v(n, 0) = r(n) / mag_r;								// try caching 1/mag_r, multiply... time diff.
			}

			for (int j = 0; j < m; j++)
			{
				sum_Av(elements, j);						// v_tilde_j+1 = A * v_j

				for (int i = 0; i <= j; i++)
				{
					h(i, j) = MPIdot(v_tilde, j+1, v, i);

					for (int n = 0; n < size; n++)
					{
						v_tilde(n, j+1) -= h(i, j) * v(n, i);
					}
				}

				h(j+1, j) = MPImag(v_tilde, j+1);

				for (int n = 0; n < size; n++)
				{
					v(n, j+1) = v_tilde(n, j+1) / h(j+1, j);					// cache 1/h[j+1][j], multiply, test?
				}


				bool stat = true;
				if (rank == 0)
				if (rotate_and_check(j+1) < tol)
				{
						cout << "mag < tol, breaking..." << endl;			// need to msg other ranks!
						stat = false;
				}
				MPI_Bcast(&stat,1,MPI_C_BOOL,0, MPI_COMM_WORLD);
				if ( !stat)
				{
					cout << "tolerence met msg received in rank " << rank << endl;
					m = j + 1;
					break;
				}
			}

			rotate();														// apply Givens rotation (Q-R algorithm)

			if (rank == 0)
				cout << "\tResidual: " << g(m) << endl;

			solve();														// solve for updates to the solution vector

			for (int n = 0; n < size; n++)									// update guess vector du
			{
				du(n) += z(n);
			}

/*			if ( abs(g(m)) < tol)											// check for convergence, restarting gmres
			{
				cout << "convergence critiera met: breaking early, iteration " << gloop << endl;
				break;
			}
*/		
			m = m_;			// restore maximum iteration size..
		}
	}

	/// Update the nodes with the residuals
	/// @param nodes NodeVar data structure, for updating results
	void update( std::vector< Node *> & nodes )
	{
		std::vector< NodeStruct > nodeList;					// the final return value, values of du correlated to places.
		std::vector< NodeStruct > nodeSendListUp;			// our request list, only used on prefetch
		std::vector< NodeStruct > nodeSendListDown;			// our request list, only used on prefetch

		for(unsigned int j=0; j<nodes.size(); j++)
		{
			int node_number = nodes[j]->number;
			if ( (node_number < min || node_number > max) )
			{
				if (!prefetchnodes)
				{
					NodeStruct dat;
					dat.index = j;
					dat.node = node_number;
					if ( node_number < min)
						nodeSendListDown.push_back( dat);				
					if ( node_number > max)
						nodeSendListUp.push_back( dat);				
				}
			}
			else
			{
				for (int k = 0; k < neqn; k++)
				{
					nodes[j]->U0(k) = nodes[j]->U(k);
					nodes[j]->U(k) += du( neqn*(node_number-delta) + k );
				}
			}
		}

		// swap data 
		if (!prefetchnodes)
		{
			MPISwapRank(nodeSendListUp, nodeReceiveListUp, rank_up, NodeStructType);
			MPISwapRank(nodeSendListDown, nodeReceiveListDown, rank_down, NodeStructType);
			prefetchnodes = true;
		}
		// now, the receive list should be in our domain, so fill it and then swap back

		for(unsigned int j = 0; j < nodeReceiveListUp.size(); j++) 
		{
			int node = (int)nodeReceiveListUp[j].node;

			for (int k = 0; k < neqn; k++)
			{
				nodeReceiveListUp[j].data[k] = du( neqn*(node-delta) + k );
			}
		}
		for(unsigned int j = 0; j < nodeReceiveListDown.size(); j++) 
		{
			int node = (int)nodeReceiveListDown[j].node;

			for (int k = 0; k < neqn; k++)
			{
				nodeReceiveListDown[j].data[k] = du( neqn*(node-delta) + k );
			}
		}

		//

		MPISwapRank(nodeReceiveListUp, nodeList, rank_up, NodeStructType);
		MPISwapRank(nodeReceiveListDown, nodeList, rank_down, NodeStructType);

		// update our nodeage from passed data
		for(unsigned int j = 0; j < nodeList.size(); j++) 
		{
			int node = (int)nodeList[j].node;
			int index = (int)nodeList[j].index;

			for (int k=0; k < neqn; k++)
			{
				nodes[index]->U0(k) = nodes[index]->U(k);
				nodes[index]->U(k) += nodeList[j].data[k];
			}
		}
	}

private:

	void prefetch(std::vector< Element *> & elements)
	{
		std::vector< DataStructAv > duSendListUp;				// our request list
		std::vector< DataStructAv > duSendListDown;				// our request list

		// try pre-calculating nodes to transmit up/down
		// rhs, Ax, Av.
		// need to store element - node combo.
		for(unsigned int j = 0; j < elements.size(); j++) 
		{
			for (int i = 0; i < nnod; i++)
			{
				int gi = elements[j]->node[i]->number; 

				if ( gi < min || gi > max )
				{
					std::vector<int> a_pair;
					a_pair += j, i, gi;
					if (gi < min)
					{
						nodesDown.push_back( a_pair);
						
						int elem =  j;
						int node =  i;
						int globnode = gi;
						DataStructAv dat;
						dat.indice = globnode;
						dat.elno = elem;
						dat.locnode = node;
						for (int k = 0; k < neqn; k++)
						{
							dat.data[k]   = 0;
						}
						duSendListDown.push_back( dat);
						
					}
					if (gi > max)
					{
						nodesUp.push_back( a_pair);

						int elem =  j;
						int node =  i;
						int globnode = gi;
						DataStructAv dat;
						dat.indice = globnode;
						dat.elno = elem;
						dat.locnode = node;
						for (int k = 0; k < neqn; k++)
						{
							dat.data[k]   = 0;
						}
						duSendListUp.push_back( dat);
					}
				}
			} 
		}
		cout << rank << "nodes up:   " << nodesUp.size() << endl;
		cout << rank << "nodes down: " << nodesDown.size() << endl;

		MPISwapRank(duSendListDown, duReceiveListDown, rank_down,DataStructAvType);
		MPISwapRank(duSendListUp,   duReceiveListUp, rank_up,DataStructAvType);
	}

		/// calculate RHS vector from element, ebe
	void compose_rhs(std::vector<Element *> & elements)
	{
		b.clear();

		std::vector< DataStruct > dataDown;
		std::vector< DataStruct > dataUp;

		for(unsigned int j = 0; j < elements.size(); j++) 
		{
			for (int i = 0; i < nnod; i++)
			{
				int gi = elements[j]->node[i]->number; 
				if ( gi >= min && gi <= max )
				{
					for(int k = 0; k < neqn; k++)
					{
						b( neqn*(gi-delta) + k ) += elements[j]->rhs(neqn*i+k);
					}
				}
			} 
		}

		// this form repeats a lot: we should either be able to 
		// (a) make a generic DataStruct in prefetch() and just fill in data.
		// (b) make a function at least to abstract this out.
		for (unsigned int i = 0; i < nodesUp.size(); i++)
		{
			int elem = nodesUp[i][0];
			int node = nodesUp[i][1];
			int globnode = nodesUp[i][2];
			DataStruct dat;
			dat.indice = globnode;

			for (int k = 0; k < neqn; k++)
			{
				dat.data[k]   = elements[elem]->rhs(neqn*node+k);
			}
			dataUp.push_back( dat);
		}
		for (unsigned int i = 0; i < nodesDown.size(); i++)
		{
			int elem = nodesDown[i][0];
			int node = nodesDown[i][1];
			int globnode = nodesDown[i][2];
			DataStruct dat;
			dat.indice = globnode;

			for (int k = 0; k < neqn; k++)
			{
				dat.data[k]   = elements[elem]->rhs(neqn*node+k);
			}
			dataDown.push_back( dat);
		}

		std::vector< DataStruct > update;

		MPISwapRank(dataDown,update,rank_down,DataStructType);
		MPISwapRank(dataUp,update,rank_up,DataStructType);

		// now take new data and add to b.
		for(unsigned int j = 0; j < update.size(); j++) 
		{
			int gi = (int)update[j].indice;
			for (int k = 0; k < neqn; k++)
			{
				double rhs = update[j].data[k];
				b( neqn*(gi-delta) + k ) += rhs;
			}
		}
	}


		/// Sum A * v, ebe
	void sum_Av(std::vector<Element *> & elements, int jj)
	{
		std::vector< DataStructAv> duList;						// the final return value, values of du correlated to places.

		// now, the receive list should be in our domain, so fill it and then swap back

		for(unsigned int j = 0; j < duReceiveListDown.size(); j++) 
		{
			int gi = (int)duReceiveListDown[j].indice;
			for (int k = 0; k < neqn; k++)
			{
				duReceiveListDown[j].data[k] = v(neqn * (gi - delta) + k, jj);
			}
		}
		for(unsigned int j = 0; j < duReceiveListUp.size(); j++) 
		{
			int gi = (int)duReceiveListUp[j].indice;
			for (int k = 0; k < neqn; k++)
			{
				duReceiveListUp[j].data[k] = v(neqn * (gi - delta) + k, jj);
			}
		}

		///

		MPISwapRank(duReceiveListDown, duList, rank_down, DataStructAvType);
		MPISwapRank(duReceiveListUp, duList, rank_up, DataStructAvType);

		///
		///

		std::vector< DataStruct > dataUp;
		std::vector< DataStruct > dataDown;

		for (unsigned int j=0; j < elements.size(); j++)
		{
			for (int y=0; y < nnod; y++)
			{
				int gy = elements[j]->node[y]->number;
				if (gy >= min && gy <= max)
				{
					for (int x=0; x < nnod; x++)
					{
						int gx = elements[j]->node[x]->number;
						if ( gx < min || gx > max )
						{
							DataStruct dat;
							dat.indice = gx;
							for (int xx = 0; xx < neqn; xx++)
							{
								dat.data[xx] = 0;		// clear dat
								for (int yy = 0; yy < neqn; yy++)
								{
									dat.data[xx] += elements[j]->R(x*neqn+xx,y*neqn+yy) * v((gy-delta)*neqn+yy, jj);
								}
							}
							if ( gx < min)
								dataDown.push_back( dat);
							if ( gx > max)
								dataUp.push_back( dat);
						}
						else
						for (int xx = 0; xx < neqn; xx++)
						{
							for (int yy = 0; yy < neqn; yy++)
							{
								v_tilde((gx-delta)*neqn+xx, jj+1) += elements[j]->R(x*neqn+xx,y*neqn+yy) * v((gy-delta)*neqn+yy, jj);
							}
						}
					}
				}
			}
		}

		for (unsigned int j = 0; j < duList.size(); j++)
		{
			int y    = (int)duList[j].locnode;
			int elno = (int)duList[j].elno;
			for (int yy = 0; yy < neqn; yy++)
			{
				double du = duList[j].data[yy];

				for (int x=0; x < nnod; x++)
				{
					int gx = elements[elno]->node[x]->number;

					if ( gx < min || gx > max )
					{
						DataStruct dat;
						dat.indice = gx;
						for (int xx = 0; xx < neqn; xx++)
						{
							dat.data[xx] = elements[elno]->R(x*neqn+xx,y*neqn+yy) * du; 
						}
						if ( gx < min)
							dataDown.push_back( dat);
						if ( gx > max)
							dataUp.push_back( dat);
					}
					else
					for (int xx = 0; xx < neqn; xx++)
					{
						v_tilde((gx-delta)*neqn+xx, jj+1) += elements[elno]->R(x*neqn+xx,y*neqn+yy) * du; 
					}
				}
			}
		}
		
		// swap data structs
	
		std::vector< DataStruct > update;

		MPISwapRank(dataUp,update,rank_up,DataStructType);
		MPISwapRank(dataDown,update,rank_down,DataStructType);

		// now take new data and add to b.
		for(unsigned int j = 0; j < update.size(); j++) 
		{
			int gx = (int)update[j].indice;
			for (int xx = 0; xx < neqn; xx++)
			{
				v_tilde( neqn*(gx-delta) + xx , jj+1) += update[j].data[xx];
			}
		}
	}

		/// sum A * x, ebe
		/// updated, may or may not work properly since first iteration, du is all zero's.
	void sum_Ax(std::vector<Element *> & elements, Tensor<double, 1> & Ax, Tensor<double, 1> & xn)
	{
		std::vector< DataStructAv > duList;						// the final return value, values of du correlated to places.

		// now, the receive list should be in our domain, so fill it and then swap back

		for(unsigned int j = 0; j < duReceiveListUp.size(); j++) 
		{
			int gi = (int)duReceiveListUp[j].indice;
			for (int k = 0; k < neqn; k++)
			{
				duReceiveListUp[j].data[k] = xn(neqn * (gi - delta) + k);
			}
		}
		for(unsigned int j = 0; j < duReceiveListDown.size(); j++) 
		{
			int gi = (int)duReceiveListDown[j].indice;
			for (int k = 0; k < neqn; k++)
			{
				duReceiveListDown[j].data[k] = xn(neqn * (gi - delta) + k);
			}
		}
		//

		MPISwapRank(duReceiveListUp, duList, rank_up, DataStructAvType);
		MPISwapRank(duReceiveListDown, duList, rank_down, DataStructAvType);

		///
		///

		std::vector< DataStruct > dataUp;
		std::vector< DataStruct > dataDown;
		for (unsigned int j=0; j < elements.size(); j++)
		{
			for (int y = 0; y < nnod; y++)
			{
				int gy = elements[j]->node[y]->number;		// swap gy / gx ordering for iteration
				if ( gy >= min && gy <= max)
				{
					for (int x = 0; x < nnod; x++)
					{
						int gx = elements[j]->node[x]->number;
						if ( gx < min || gx > max )
						{
							DataStruct dat;
							dat.indice = gx;

							for (int xx = 0; xx < neqn; xx++)
							{
								dat.data[xx] = 0;		// clear dat
								for (int yy = 0; yy < neqn; yy++)
								{
									dat.data[xx] += elements[j]->R(x*neqn+xx,y*neqn+yy) * xn((gy-delta)*neqn+yy);
								}
							}
							if ( gx < min)
								dataDown.push_back( dat);
							else if ( gx > max)
								dataUp.push_back( dat);
						}
						else
						for (int xx = 0; xx < neqn; xx++)
						{
							for (int yy = 0; yy < neqn; yy++)
							{
								Ax((gx-delta)*neqn+xx) += elements[j]->R(x*neqn+xx,y*neqn+yy) * xn((gy-delta)*neqn+yy);
							}
						}
					}
				}
			}
		}

		for (unsigned int j = 0; j < duList.size(); j++)
		{
			int y    = (int)duList[j].locnode;
			int elno = (int)duList[j].elno;

			for (int yy = 0; yy < neqn; yy++)
			{
				double du = duList[j].data[yy];

				for (int x=0; x < nnod; x++)
				{
					int gx = elements[elno]->node[x]->number; 
					if ( gx < min || gx > max )
					{
						DataStruct dat;
						dat.indice = gx;
						for (int xx = 0; xx < neqn; xx++)
						{
							dat.data[xx] = elements[elno]->R(x*neqn+xx,y*neqn+yy) * du; 
						}
						if ( gx < min)
							dataDown.push_back( dat);
						if ( gx > max)
							dataUp.push_back( dat);
					}
					else
					for (int xx = 0; xx < neqn; xx++)
					{
							Ax((gx-delta)*neqn+xx) += elements[elno]->R(x*neqn+xx,y*neqn+yy) * du;
					}
				}
			}
		}

		// swap data structs	
		std::vector< DataStruct > update;

		MPISwapRank(dataUp,update, rank_up,DataStructType);
		MPISwapRank(dataDown,update, rank_down,DataStructType);

		// now take new data and add to b.
		for(unsigned int j = 0; j < update.size(); j++)			/// VALID
		{
			int gx = (int)update[j].indice;
			for (int xx = 0; xx < neqn; xx++)
			{
				double rhs = update[j].data[xx];
				Ax( neqn*(gx-delta) + xx ) += rhs;
			}
		}
	}

		// Turn by turn rotation to check convergence
	double rotate_and_check(int n)
	{
		Tensor<double, 2> hh(n+1, n);
		Tensor<double, 1> gg(n+1);
		gg(0) = g(0);

		for (int i = 0; i < n+1; i++)
		{
			for (int j = 0; j < n; j++)
			{
				hh(i, j) = h(i, j);
			}
		}

		for (int i = 0; i < n; i++)
		{
			double mag = sqrt( pow(h(i, i), 2) + pow(h(i+1, i), 2));			// minimal effect, but can calc 1/mag ... 
			double c = h(i, i)   / mag; 
			double s = h(i+1, i) / mag; 

			for(int j=0; j < n; j++)
			{
				double h0 = hh(i, j);
				double h1 = hh(i+1, j);
				hh(i, j)   =  c * h0 + s * h1;
				hh(i+1, j) = -s * h0 + c * h1;
			}

			gg(i+1) = -s * gg(i);
			gg(i)   =  c * gg(i);
		}

		cout << "residual: " << gg(n) << endl;
		return abs(gg(n));
	}

		/// Apply givens rotation to Hessenberg matrix
	void rotate()
	{
		for (int i = 0; i < m; i++)
		{
			double mag = sqrt( pow(h(i, i), 2) + pow(h(i+1, i), 2));			// minimal effect, but can calc 1/mag ... 
			double c = h(i, i)   / mag; 
			double s = h(i+1, i) / mag; 

			for(int j=0; j < m; j++)
			{
				double h0 = h(i, j);
				double h1 = h(i+1, j);
				h(i, j)   =  c * h0 + s * h1;
				h(i+1, j) = -s * h0 + c * h1;
			}

			g(i+1) = -s * g(i);
			g(i)   =  c * g(i);
		}
	}

		/// Solve the Hessenberg system of equations by back substitution,
		/// then solve for the residual updates to du, (z)
	void solve()
	{
			// solve the Hessenberg system of equations by back-substitution
		for(int j= m-1; j>=0; j--)
		{
			y(j) = g(j);

			for(int i = m-1; i >= 0; i--)
			{
				if(i!=j)
				{
					y(j) -= h(i, j) * y(i);
				}
			}
			y(j) /= h(j, j);
		}

		for (int j=0; j<size; j++)
		{
			for (int i=0; i < m; i++)
			{
				z(j) += y(i) * v(j, i);
			}
		}
	}

		/// Could probably be optomized ...  using the mutual swapping etc. available in MPI
	template< class T>
	void MPISwapRank(std::vector<T> & data, std::vector<T> & update, int other_rank, MPI_Datatype & type)
	{
		int size = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		if (other_rank == -1 || other_rank == size)
			return;

		int tag = rank * other_rank;

		//// SEND BLOCK ////
		int num_msg_send = data.size();
		MPI_Send( &num_msg_send, 1, MPI_INT, other_rank, tag, MPI_COMM_WORLD);

		//// RECV BLOCK ////
		T dat;
		int num_msg_recv = 0;
		MPI_Recv( &num_msg_recv, 1, MPI_INT, other_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int j = 0; j < num_msg_recv; j++)
		{
			MPI_Recv( &dat, 1, type, other_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			update.push_back( dat);
		}
		//// RECV BLOCK ////

		for (int j = 0; j < num_msg_send; j++)
		{
			MPI_Send( &data[j], 1, type, other_rank, tag, MPI_COMM_WORLD);
		}
		//// SEND BLOCK ////

	}


		/// clear all data; call at the beginning of each GMRES loop
		/// this does NOT clear b or du, both of which persist through the whole iteration
		/// b is cleared in the function which updates b.
		/// du is cleared at the beginning of the gmres iterate() function
	void clear()
	{
		r.clear();							// initial residual;
		v.clear();							// 
		v_tilde.clear();					//
		h.clear();							// Hessenberg matrix, prior to Givens rotation
		g.clear();							// Ax+b to Hessenberg matrix
		y.clear();
		z.clear();							// final correction
		Ax.clear();							// A * (x_0 + x)
	}

	int restart;						/// maximum number of restarts
	int m;								/// maximum number of gmres iterations per restart (varies with flexible gmres)
	int m_;								/// maximum number of gmres iterations per restart
	int neqn;							/// number of equations
	int nnod;							/// number of nodes
	int ndim;							/// number of dimensions
	int size;							/// size of dims (nnod in domain * neqn)

	int rank;							/// our MPI rank
	int numthreads;						/// number of MPI threads (size)
	int delta;							/// thread 2(+) offset.
	int min;							/// min node, filter
	int max;							/// max node, filter
	int rank_down;						/// rank down
	int rank_up;						/// rank up

	MPI_Datatype DataStructAvType;		/// data struct Av/Ax mult type
	MPI_Datatype DataStructType;		/// data struct b type
	MPI_Datatype NodeStructType;		/// node data xfer

	bool do_once;
	bool prefetchnodes;
	std::vector< std::vector<int> > nodesUp;
	std::vector< std::vector<int> > nodesDown;
	std::vector< NodeStruct > nodeReceiveListUp;		// our list to fulfill
	std::vector< NodeStruct > nodeReceiveListDown;		// our list to fulfill
	std::vector< DataStructAv > duReceiveListUp;			// our list to fulfill
	std::vector< DataStructAv > duReceiveListDown;			// our list to fulfill

	Tensor<double, 1> r;						/// initial residual;
	Tensor<double, 2> v;						/// 
	Tensor<double, 2> v_tilde;					///
	Tensor<double, 2> h;						/// Hessenberg matrix, prior to Givens rotation
	Tensor<double, 1> g;						/// Ax+b to Hessenberg matrix
	Tensor<double, 1> y;						/// solution to the Hessenberg system of equations
	Tensor<double, 1> z;						/// final correction
	Tensor<double, 1> Ax;						/// A * (x_0 + x)
	Tensor<double, 1> du;						/// delta-u, initial guess vector
	Tensor<double, 1> b;						/// rhs of the system of equations
};

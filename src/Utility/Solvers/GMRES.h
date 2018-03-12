#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "GMRESbase.h"
#include "../Tensor.h"

	/// the Generalized Minimial RESidual method for solving
	/// the matrix problem Ax = b; this is accomplished using
	/// an element-by-element (EBE) assembly of A and b.

struct GMRES : public GMRESbase
{
	/// Constructor, intializes dimensions of arrays
	GMRES(int r, int m, int size, int neqn, int nnod, int ndim)
		: size(size),
		  restart(r), 
		  m(m),
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
		  ndim(ndim)
	{
		cout << "GMRES unparallelized" << endl;
		ofstream fout("gmres.csv",ios::out);
		fout.close();
	}

	/// Solves the system of equations, given the element correlation and node data
	/// @param elements element data structure, provides input data
	void iterate( std::vector< Element *> & elements )
	{
		cout << "GMRES::Iterate" << endl;
		ofstream fout("gmres.csv",ios::app);

		du.clear();
		
		compose_rhs(elements);				// clear, then compose b

		double tol = mag(b) * 1.0E-9;		// calculate tolerence for convergence (to quit gmres restarts early)

		fout << "tol, " << tol << ",res:";

		cout << "tolerence: " << tol << endl;
		cout << "residual: ";

		for (int gloop = 0; gloop < restart; gloop++)
		{
			clear();						// clear nonpersistant data -> PROBLEM WITH THIS CALL

			sum_Ax(elements,Ax,du);
			// sum initial residual
			for (unsigned int n = 0; n < r.imax(); n++)
			{
				r(n) = b(n) - Ax(n);
			}

			// calculate RMS of residual
			double mag_r = mag(r);

			g(0) = mag_r;			// intialize Hessenberg matrix solution vector

			for (int n = 0; n < size; n++)
			{
				v(n, 0) = r(n) / mag_r;								// try caching 1/mag_r, multiply... time diff.
			}

			for (int j=0; j<m; j++)
			{
				sum_Av(elements, j);						// v_tilde_j+1 = A * v_j

				for (int i=0; i<=j; i++)
				{
					h(i, j) = dot(v_tilde, j+1, v, i);

					for (int n = 0; n < size; n++)
					{
						v_tilde(n, j+1) -= h(i, j) * v(n, i);
					}
				}

				h(j+1, j) = mag(v_tilde, j+1);

				for (int n = 0; n < size; n++)
				{
					v(n, j+1) = v_tilde(n, j+1) / h(j+1, j);					// cache 1/h[j+1][j], multiply, test?
				}

				/*
				rotate(j);
				if (abs(g(j+1)) < tol)
				{
					cout << "tolerence met, i = " << j << ", breaking, " << g(j+1) << endl;
					break;
				}
					// need to verify the rotations are OK, then solve just the subset.
				*/
			}

			// apply Givens rotation (Q-R algorithm)
			rotate();

			fout << ", " << g(m);

			cout << g(m) << ", ";		// outputing original residual?

			// solve for updates to the solution vector
			solve();

			// update guess vector delta-U
			for (int n = 0; n < size; n++)
			{
				du(n) += z(n);
			}

			// check for convergence
			if ( abs(g(m)) < tol)
			{
				cout << "convergence critiera met: breaking early, iteration " << gloop << endl;
				break;
			}
		}

		fout << endl;
		cout << endl;
		fout.close();
	}

	/// Update the nodes with the residuals
	/// @param nodes NodeVar data structure, for updating results
	void update( std::vector< Node *> & nodes )
	{
		for(unsigned int j=0; j<nodes.size(); j++)
		{
			for (unsigned int k=0; k<nodes[j]->U.imax(); k++)
			{
				nodes[j]->U0(k) = nodes[j]->U(k); 
				nodes[j]->U(k) += du( neqn*j + k );
			}
		}
	}

private:

	//
	//	In one thread, putting back in the "right order" the following three functions'
	//	orderings increases GMRES speed by 17%. This makes sense as we aren't jumping 
	//	around in memory so much. If the end goal is MPI parallelization, then we 
	//	shouldn't care at this point about the OpenMP parallelization ordering.
	//

		/// calculate RHS vector from element, ebe
	void compose_rhs(std::vector<Element *> & elements)
	{
		// b = 0;
		b.clear();

		for(unsigned int j = 0; j < elements.size(); j++) 
		{
			for (int i = 0; i < nnod; i++)
			{
				int gi = elements[j]->node[i]->number; 
				for(int k = 0; k < neqn; k++)
				{
					b( neqn*gi + k ) += elements[j]->rhs(neqn*i+k);
				}
			} 
		}
	}
		/// Sum A * v, ebe
	void sum_Av(std::vector<Element *> & elements, int jj)
	{
		for (unsigned int j=0; j < elements.size(); j++)
		{
			for (int x=0; x < nnod; x++)
			{
				int gx = elements[j]->node[x]->number;
				for (int y=0; y < nnod; y++)
				{
					int gy = elements[j]->node[y]->number;

					for (int xx = 0; xx < neqn; xx++)
					{
						for (int yy = 0; yy < neqn; yy++)
						{
							v_tilde(gx*neqn+xx, jj+1) += elements[j]->R(x*neqn+xx,y*neqn+yy) * v(gy*neqn+yy, jj);
						}
					}
				}
			}
		}
	}

		/// sum A * x, ebe
	void sum_Ax(std::vector<Element *> & elements, Tensor<double, 1> & Ax, Tensor<double, 1> & xn)
	{
		for (unsigned int j=0; j < elements.size(); j++)
		{
			for (int x = 0; x < nnod; x++)
			{
				int gx = elements[j]->node[x]->number;
				for (int y = 0; y < nnod; y++)
				{
					int gy = elements[j]->node[y]->number;
					for (int xx = 0; xx < neqn; xx++)
					{
						for (int yy = 0; yy < neqn; yy++)
						{
							Ax(gx*neqn+xx) += elements[j]->R(x*neqn+xx,y*neqn+yy) * xn(gy*neqn+yy);
						}
					}
				}
			}
		}
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

		/// Apply givens rotation to Hessenberg matrix
	void rotate(int i)
	{
//		for (int i = 0; i < m; i++)
//		{
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
//		}
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
	int m;								/// maximum number of gmres iterations per restart
	int neqn;							/// number of equations
	int nnod;							/// number of nodes
	int ndim;							/// number of dimensions
	int size;							/// size of dims (nnod in domain * neqn)

	Tensor<double, 1> r;							/// initial residual;
	Tensor<double, 2> v;							/// 
	Tensor<double, 2> v_tilde;					///
	Tensor<double, 2> h;							/// Hessenberg matrix, prior to Givens rotation
	Tensor<double, 1> g;							/// Ax+b to Hessenberg matrix
	Tensor<double, 1> y;							/// solution to the Hessenberg system of equations
	Tensor<double, 1> z;							/// final correction
	Tensor<double, 1> Ax;						/// A * (x_0 + x)
	Tensor<double, 1> du;						/// delta-u, initial guess vector
	Tensor<double, 1> b;							/// rhs of the system of equations
};

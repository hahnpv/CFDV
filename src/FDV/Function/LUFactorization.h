#pragma once

#include <iostream>
using namespace std;

	/// Performs L-U Factorization for LU-Preconditioning
	/// does not actually, at this time, perform preconditioning
	/// maybe we can have it do the first hit, but within GMRES, 
	/// or another solver, we need to figure out the best way how
	/// before implementing.

	/// Has NOT been verified.
template<class T> struct LUFactorization : public unary_function<T, void>
{
	LUFactorization(int nnod, int neqn)
		: dim(nnod*neqn)
	{
	}

	void operator() (T e) 
	{
		// Clear L and U since ClearElement does not do this by default
		// (preconditioning is currently an option, not de-facto)

		e->L.clear();
		e->U.clear();

		// using the Doolittle method, pg 896 of Kreyszig, Advanced
		// Engineering Mathematics, L=m_jk and U=u_jk

		for (int k = 0; k < dim; k++)
		{
			e->U(0, k) = e->R(0, k);
		}

		for (int j = 1; j < dim; j++)
		{
			e->L(j, 0) = e->R(j, 0) / e->U(0, 0);
		}

		for (int j = 1; j < dim; j++)
		{
			for (int k = j; k < dim; k++)
			{
				e->U(j, k) = e->R(j, k);

				for (int s = 0; s < (j - 1); s++)
				{
					e->U(j, k) -= e->L(j, s) * e->U(s, k);
				}
			}
		}

		for (int k = 1; k < dim; k++)
		{
			for (int j = (k + 1); j < dim; j++)
			{
				double mu = 0;
				for (int s = 0; s < (k - 1); s++)
				{
					mu += e->L(j, s) * e->U(s, k);
				}
				e->L(j, k) = (e->R(j, k) - mu) / e->U(k, k);
			}
		}
	}
	const int dim;
};


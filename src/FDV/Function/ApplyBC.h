#pragma once

	/// Old version (as of 01/2009) was not working properly in 3D. Removing it
	/// from consideration made the code work. Re-working ApplyBC to function properly.
	/// works allright with a single process, breaks in multiple threads?
	/// residual breaks halfway through iteration 2?
template<class T> struct ApplyBC : public unary_function<T, void>
{
  ApplyBC(int neqn, int nnod, int ndim, int nbnod)
	  : neqn(neqn),
		nnod(nnod),
		ndim(ndim),
		nbnod(nbnod)
	{
/*		if (ndim == 2)
			prescribe = 0.5;
		else if (ndim == 3)
			prescribe = 0.25;
*/	}

	void operator() (T e) 
	{
		Thermo & thermo = Thermo::Instance();
		double Cv = thermo.Cv;
		double Twall = thermo.Twall;
/*
		for (int f = 0; f < e->face.size(); f++)
		{
			for (int i = 0; i < nbnod; i++)
			{
				int ni = e->face[f]->n(i);
				for (int j = 0; j < neqn; j++)		// skip E for now ... E not the problem... rhs is.
				{
					if (e->node[ ni ]->dirichlet(j))
					{
						for (int ii = 0; ii < nnod*neqn; ii++)
						{
//							e->R(ii, i*neqn + j) = 0;
							e->R(ni*neqn + j, ii) = 0;
						}
						for (int ii = 0; ii < nbnod; ii++)
						{
							for (int jj = 0; jj < neqn; jj++)
							{
								e->R( ni * neqn + jj, ni * neqn + j) = 0;
							}
						}

						e->R(ni*neqn + j, ni*neqn + j) = e->node[ni]->coeff;
						e->rhs(ni*neqn + j) = e->node[ni]->coeff * e->node[ni]->U0(j);

						if (j == neqn-1 && !e->node[ni]->dirichlet(0) )			// verify this term is ok. (flat plate)
						{
							e->rhs(ni*neqn + j) = 0.;
	//						e->R(i*neqn+j, i*neqn) = -Cv * Twall * e->node[i]->coeff;
						}
					}	
				}
			}
		}
*/

		for (int i = 0; i < nnod; i++)
		{
			for (int j = 0; j < neqn; j++)		// skip E for now ... E not the problem... rhs is.
			{
				if (e->node[i]->dirichlet(j))
				{
					for (int ii = 0; ii < nnod*neqn; ii++)
					{
						e->R(ii, i*neqn + j) = 0;
						e->R(i*neqn + j, ii) = 0;
					}

					e->R(i*neqn + j, i*neqn + j) = e->node[i]->coeff;
					e->rhs(i*neqn + j) = e->node[i]->coeff * e->node[i]->U0(j);

					if (j == neqn-1 && !e->node[i]->dirichlet(0) )			// verify this term is ok. (flat plate)
					{
						e->rhs(i*neqn + j) = 0.;
//						e->R(i*neqn+j, i*neqn) = -Cv * Twall * e->node[i]->coeff;
					}
				}	
			}
		}

	}

	int neqn;
	int nnod;
	int ndim;
	int nbnod;

//	double prescribe;
};


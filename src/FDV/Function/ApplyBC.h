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

		// Iterate over faces
		for (int f = 0; f < e->face.size(); f++)
		{
			// FIXME this never hits!
			for (int inode = 0; inode < nbnod; inode++)									// face node iterator 1
			{
				int ia = e->face[f]->n(inode);
				for (int ieqn = 0; ieqn < neqn; ieqn++)									// face node equation iterator 1
				{
					int irow = neqn * ia + ieqn;
					for (int jnode = 0; jnode < nbnod; jnode++)							// face node iterator 2
					{
						int ja = e->face[f]->n(jnode);
						for (int jeqn = 0; jeqn < neqn; jeqn++)							// face node equation iterator 2
						{
							int icol = neqn * ja + jeqn;

							// dirichlet. Note: this is a face property not a nodal one. FIXME
							if (e->node[ia]->dirichlet(ieqn))
							{
		/*						// If Wall, [this is not working] FIXME doesnt work / matrix looks f'd
								if (e->face[f]->bc == 1 && ieqn==neqn-1	)
								{
									int jcol = icol - (ndim + 1);							// points to rho when you are on energy term
									if (irow != icol)
									{
										e->R(irow, jcol) += e->R(irow, icol) * Cv * Twall;
										e->R(irow, icol) = 0.;
									}
									else
									{
										for (int k = 0; k < neqn*nnod; k++)
										{
											e->R(irow, k) = 0.;
										}
										e->R(irow, jcol) = -Cv * Twall * 0.5;
										e->R(irow, icol) = 0.5;
										e->rhs(irow) = 0.;
									}
								}
	*/							
								// If non-wall 
//								if (e->face[f]->bc == -1)			// this appears to be working! inflow
								{
									e->R(irow, icol) = 0.;
									e->rhs(irow) = 0.;
									if (irow == icol)
									{
										for (int k = 0; k < neqn*nnod; k++)
										{
											e->R(irow, k) = 0.;
										}
										e->R(irow, icol) = 0.5;
									}
								}
							}
						}
					}
				}
			}
			/*
			if (e->face[f]->bc == 1)
			{
				cout << e->R << endl << endl;
				cout << e->rhs << endl << endl;
				cout << endl;
			}
			*/
		}
	}

	int neqn;
	int nnod;
	int ndim;
	int nbnod;

//	double prescribe;
};


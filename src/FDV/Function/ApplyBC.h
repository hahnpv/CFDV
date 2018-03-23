#pragma once

	/// Being reworked
template<class T> struct ApplyBC : public unary_function<T, void>
{
  ApplyBC(int neqn, int nnod, int ndim, int nbnod)
	  : neqn(neqn),
		nnod(nnod),
		ndim(ndim),
		nbnod(nbnod),
		first_pass(true)
	{}

	void operator() (T e) 
	{
		Thermo & thermo = Thermo::Instance();
		double Cv = thermo.Cv;
		double Twall = thermo.Twall;

		// Iterate over faces
		for (int f = 0; f < e->face.size(); f++)
		{
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
							// but carefully because edge of plates need handling!
							if (e->node[ja]->dirichlet(jeqn))
							{
								// If Wall
								if (e->face[f]->bc == 1 && jeqn==neqn-1	&& !thermo.adiabatic)
								{/*
									// could problem be edges of plate where we null out local but not global?
									// ACTUALLY should be fine because node says dirichlet even if bc does not
									int jcol = icol - (ndim + 1);						// points to rho when you are on energy term
									if (irow != icol)
									{
										e->R(jcol, irow) += e->R(icol, irow) * Cv * Twall; // jcol,irow is E not rho in mat?
										e->R(icol, irow) = 0.;
									}
									else
									{
										e->rhs(icol) = 0.;
										for (int k = 0; k < neqn*nnod; k++)
										{
											e->R(icol, k) = 0.;
										}
//										cout << jcol << " " << irow << " " << icol << endl; 4,7,7, 0,3,3
										e->R(jcol, irow) = -Cv * Twall*0.5;			// two contributions not sure if this is right
										e->R(icol, irow) = 0.5;					// why is 0.5 needed?
									}
									*/
								}
								else // this works with isothermal logic above is ignored
								{
									e->R(icol, irow) = 0.; // FIXME make sure you understand col/row tensor
									e->rhs(icol) = 0.;
									if (irow == icol)
									{
										for (int k = 0; k < neqn*nnod; k++)
										{
											e->R(icol, k) = 0.;
										}
										e->R(icol, irow) = 0.5;
									}
								}


								// could problem be edges of plate where we null out local but not global?
								// ACTUALLY should be fine because node says dirichlet even if bc does not
								/*
								int jcol = icol - (ndim + 1);						// points to rho when you are on energy term
								if (irow != icol)
								{
									e->R(jcol, irow) += e->R(icol, irow) * Cv * Twall; // jcol,irow is E not rho in mat?
									e->R(icol, irow) = 0.;
								}
								else
								{

									for (int k = 0; k < neqn*nnod; k++)
									{
										e->R(icol, k) = 0.;
									}
									e->R(jcol, irow) = -Cv * Twall;
									e->R(icol, irow) = 0.5;
									e->rhs(icol) = 0.;

								}
								*/
								// This works for fixed condx (rho, u, v, rho+E)
								// This does not work for isothermal wall (rho+E)
								// this works with e->node[ia]->dirichlet(ieqn), using i as iterator
								// above inversion works with j as iterator
/*									e->R(irow, icol) = 0.;
									e->rhs(irow) = 0.;
 									if (irow == icol)
									{
										for (int k = 0; k < neqn*nnod; k++)
										{
											e->R(irow, k) = 0.;
										}
										e->R(irow, icol) = 0.5;
									}
*/
							}
						}
					}
				}
			}
/*
			if (e->face[f]->bc == 1 && first_pass)
			{
				cout.precision(4);
				cout << e->R << endl << endl;
				cout << e->rhs << endl << endl;
				cout << endl;
//				cin.get();
				first_pass = false;
			}
*/
		}
	}

	int neqn;
	int nnod;
	int ndim;
	int nbnod;

	bool first_pass;		// debug
};


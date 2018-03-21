#pragma once

	/// Enforce Dirichlet boundary conditions by reapplying
	/// prior value of the conservation variable, U0
template<class N> struct EnforceBC : public unary_function<N, void>
{
	EnforceBC(int ndim)
	{
		this->ndim = ndim;
//		cout << "EnforceBC no longer filters out Node[]->bc of 0 or >100, verify ..." << endl;
		// this works, for very large cases might waste a few seconds ... think if thre is something better.
	}
	void operator() (N n) 
	{
	//	if (e->bc != 0 && e->bc < 100)				// Disregard non-boundaries and flagged conditions (not true boundaries)
	//	{
//		if (e->bc != 1)
//		{
			if (n->dirichlet(0))
			{
				n->U(0) = n->U0(0);
			}

			for (int i=1; i < ndim+1; i++)			// rho, v[]
			{
				if (n->dirichlet(i))
				{
//					n->U(i) = n->U0(i);
					n->U(i)  = ((n->U0(i)) / (n->U0(0))) * n->U(0);
				}
			}
			if (n->dirichlet(ndim+1))				// E
			{
				n->U(ndim+1)  = ((n->U0(ndim+1)) / (n->U0(0))) * n->U(0);		
			}
//		}
	//	}
/*		if (e->bc == 1)
		{
			int neqn = e->U.imax;
			BoundaryConditionInflow inflow(neqn);
			inflow.ApplyDirichlet( *e);
		}
*/	}

	int ndim;
};


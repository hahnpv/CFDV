#pragma once

#include <iostream>
using namespace std;

	/// Calculates an updated timestep using the CFL criteria
template<class T> struct CFL : public unary_function<T, void>
{
	CFL(std::vector<T> & elements, double & cfl, double & delta_t)
		: cfl(cfl)
	{
		dx_dv = 1000;	// high initial guess
 		for(unsigned int i=0; i<elements.size(); i++)
		{
			(*this)(elements[i]);				
		}

		MPI::COMM_WORLD.Allreduce(&dt, &delta_t, 1, MPI_DOUBLE, MPI_MIN);
	}

	void operator() (T e) 
	{
		Thermo & thermo = Thermo::Instance();

		double dxele = 1000;

		for (unsigned int i = 0; i < e->node.size(); i++)
		{
			double Temp = 0;
		
			if(e->node[i]->T > 10E-10)
			{
				Temp = e->node[i]->T;
			}

			double v = mag(e->node[i]->v);		

			double a = sqrt(1.4*thermo.cgas*Temp);

			if (v+a > 10E-10)						// avoid division by zero
			{
				double dxnode = e->cleng / ( v + a);

				if (dxnode < dxele)
				{
					dxele = dxnode;
				}
			}
		}

		if ( dxele < dx_dv)
		{
			dx_dv = dxele;
			dt = cfl * dx_dv;
		}
	}

  double cfl;
  double dx_dv;
  double dt;
};

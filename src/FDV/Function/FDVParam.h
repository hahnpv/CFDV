#pragma once

// Calculate the FDV parameters

template<class E> struct FDVParam : public unary_function<E, void>
{
	FDVParam()
	{
		eta = 0.2;								// 0.05 - 0.2
	}

	void operator() (E e, double & s1, double & s2, double & s3, double & s4) 
	{
		Thermo & thermo = Thermo::Instance();
	
		double r = 0;
			
		double M_min = 0;
		double M_max = 0;
		double Re_min = 0;
		double Re_max = 0;

		for (unsigned int beta=0; beta < e->node.size(); beta++)
		{
			double V = mag(e->node[beta]->v);
			double T = e->node[beta]->T;

			if ( T < 1.0E-12)
			{
				T = 0;
			}

			double a = sqrt( thermo.gamma * thermo.cgas * T );
			double M = 0;
			if (a <1.0E-12)			// mean zero or negative temperature.
			{
			//	M = thermo.cmach;
				M = 0;				//
			}
			else
			{
				M = V / a;
			}

			double mu = thermo.calc_mu( T);

			double Re = (e->node[beta]->rho * V) / mu;

			if ( beta == 0)			// on the first node set max and min to M to initialize.
			{
				M_max = M;
				M_min = M;

				Re_max = Re;
				Re_min = Re;
			}

			M_max = max(M_max,M);
			M_min = min(M_min,M);

			Re_max = max(Re_max,Re);
			Re_min = min(Re_min,Re);
		}

		// s1,s2 
		double eps = 0.0000001;
		r = sqrt( pow(M_max, 2) - pow(M_min, 2) ) / (M_min+eps); // FIXME

		if ( r > 0.01)
		{
//			e->s1 = min(r, 1.0);
			s1 = min(r, 1.0);
		}
		else if ( r < 0.01 && M_min >= 10E-10)
		{
//			e->s1 = 0.0;
			s1 = 0;
		}
		else if (M_min < 10E-10)			// catch damn-near-zero's (ie, it's zero, but roundoff in squaring/rooting)
		{
//			e->s1 = 1.0;
			s1 = 1.0;
		}

//		e->s2 = 0.5 * ( 1 + pow(e->s1, eta));
		s2 = 0.5 * ( 1 + pow(s1, eta));

		// s3,s4
		r = sqrt( pow(Re_max, 2) - pow(Re_min, 2) ) / (Re_min+eps); // FIXME

		if ( r > 0.01)
		{
//			e->s3 = min(r, 1.0);
			s3 =  min(r, 1.0);
		}
		else if (r < 0.01 && Re_min >= 10E-10)
		{
//			e->s3 = 0.0;
			s3 = 0.0;
		}
		else if (Re_min <= 10E-10)
		{
//			e->s3 = 1.0;
			s3 = 1.0;
		}

//		e->s4 = 0.5 * ( 1 + pow(e->s3, eta));
		s4 =  0.5 * ( 1 + pow(s3, eta));
	}

	double eta;
};


#pragma once
#include "../Utility/Singleton.h"
#include <iostream>
using namespace std;

	// Utilizes a templated singleton
class Thermo : public Singleton<Thermo>
{
public:

	void setGamma( double g)
	{
		gamma = g;
		gamm1 = (gamma - 1.0);
		gm1d2 = ((gamma - 1.0)/2.0);
		gamm3 = (gamma - 3.0);
		gm3d2 = ((gamma - 3.0)/2.0);
	}

	double gamma;
	double gamm1;
	double gm1d2;
	double gamm3;
	double gm3d2;

	// these values can fluxuate, but constant within one run
	double Pr;									// prescribed Prandtl number
	double csuth;								// nondim. Sutherland's law for mu (using Tinf from file)
	double Cv;									// nondim. 
	double cgas;								// nondim gas constant
	double cmach;								// mach number
	double creyn;								// reynolds number
	double Twall;								// wall temperature if isothermal

	bool adiabatic;								// adiabatic wall condition (true/false)

	double calc_mu( double T)
	{
		if (T < 1.0E-12)
		{
			return 1.0;
		}
		else
		{
		//  OPTOMIZATION: T^1.5 = T*T^0.5 = T*sqrt(T)
		//	Under gcc for some reason pow(T,1.5)'s execution is nonlinear, ie,
		//	it takes a very long time for awhile but later in execution it is quick
		//	sqrt(T) is faster than pow(T) and does not exhibit this tendancy.
		//	Problem was not noted under MSVC.
		//	Results will be identical.
		//	return ( (1+csuth)/(T + csuth))*pow(T,1.5);
			return ( (1+csuth)/(T + csuth))*T*sqrt(T);
		}
	}

	double calc_k(double mu)
	{
		return mu / ( (gamma-1.0)*(cmach*cmach)*Pr );
	}

protected: 
    Thermo()									// default constructor 
	{}
private:
    friend class Singleton<Thermo>;
//	Singleton(); // ctor hidden
	Thermo& operator=(Thermo const&); // assign op. hidden
//	~Thermo(); // dtor hidden
	Thermo(Thermo const&); // copy ctor hidden -> need another method to enable copy constructor
};

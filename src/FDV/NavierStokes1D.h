#pragma once

#include "Flow.h"
#include "Thermo.h"
#include <cmath>
#include "NavierStokes.h" 
struct NavierStokes1D : public NavierStokes
{
	// 2D Convective Flux
	void calc_F(Tensor<double, 1> & U, double p, Tensor<double, 2> & matrix)
	{
		matrix.clear();

		double rho = U(0);
		double l = U(1);
		double m = U(2);
		double e = U(3);

		matrix(0, 0) = l;
		matrix(0, 1) = p + pow(l,2)/rho;
		matrix(0, 2) = (p + e)*(l/rho);
	}

	// 2D Viscous Flux
	void calc_G(Tau & tau, Tensor<double, 1> & v, Tensor<double, 1> & dT, double & HK, double & Re, Tensor<double, 2> & matrix)
	{	
		matrix.clear();

		double reinv = 1.0 / Re;

		matrix(0, 0) =0.0;
		matrix(0, 1) =-tau(0, 0)*reinv;
		matrix(0, 3) =(-tau(0, 0)*v(0) - HK*dT(0))*reinv;
	}

	/// 2D Source Term, not incorportated
	void calc_B(Node * n, Tensor<double,1> & matrix)
	{
		matrix.clear();
	};


		// 2D Flux Jacobians

	// dF/dU - verified w/GWH
	void calc_a(Flow & f, Tensor<double,3> & matrix)
	{
		Thermo & thermo = Thermo::Instance();
		matrix.clear();

      double u  = f.v(0);
      double uu = u*u;
	  double E = f.E;

	  matrix(0, 0, 02) = 0.0;
      matrix(0, 0, 1) = 1.0;
	  matrix(0, 0, 2) = 0.0;

	  matrix(0, 1, 0) = thermo.gm3d2*uu;
	  matrix(0, 1, 1) =-thermo.gamm3*u;
	  matrix(0, 1, 2) = thermo.gamm1;

	  matrix(0, 2, 0) =-thermo.gamma*E*u+thermo.gamm1*u*(uu);
	  matrix(0, 2, 1) = thermo.gamma*E-thermo.gm1d2*(3.*uu);
      matrix(0, 2, 2) = thermo.gamma*u;
	}

	// dG/dU
	void calc_b(Flow & flow, Tensor<double, 3> & matrix)
	{
		Thermo & thermo = Thermo::Instance();

		matrix.clear();

		double reinv = 1.0/thermo.creyn;
		double mu = flow.mu;
		Tensor<double, 2> & tau = flow.tau;
		double E = flow.E;

      double lambda=-2.*mu/3.;
      double mu_R=2.*mu+lambda;
      double DENOM=1./pow(flow.U(0),2);
	  double TERM1=flow.k*DENOM/thermo.Cv;

      double rho = flow.U(0);
      double u = flow.v(0);
      double uu = u*u;

	  double drdx = flow.dU(0, 0);
      double dldx = flow.dU(1, 0);
      double dedx = flow.dU(2, 0);

		// check ordering of tau, etc.

	  matrix(0, 1, 0) =(DENOM*(mu_R*dldx-2.*(mu_R*u*drdx)))*reinv;
      matrix(0, 1, 1) =mu_R*drdx*DENOM*reinv;

	  // may be issue with a few terms here, third term might need a minus sign?
	  matrix(0, 2, 0) =u*matrix(0, 1, 0)+((u*tau(0,0))/rho -TERM1*(-dedx+(2.*E-3.*uu)*drdx+2.*u*dldx))*reinv;
      matrix(0, 2, 1) =u*matrix(0, 1, 1)+(-tau(0,0)/rho-TERM1*(2.*u*drdx-dldx))*reinv;
      matrix(0, 2, 3) =TERM1*drdx*reinv;
	}
//	Thermo thermo;


	// dG/dU,j		
	void calc_c(Flow & flow, Tensor<double, 4> & matrix)
	{
		Thermo & thermo = Thermo::Instance();
		matrix.clear();

		double reinv = 1.0/thermo.creyn;
		double mu = flow.mu;

      double lambda=-2.*mu/3.;
      double mu_R=2.*mu+lambda;

      double rho      = flow.U(0);
      double denom  = 1./pow(rho,2);
      double denom2 = denom/rho;

      double l     = flow.U(1);
      double re      = flow.U(2);
      double ll    = l * l;
	  double hkocv   = flow.k/thermo.Cv;

	  // good to here

      matrix(0, 0, 1, 0) =  mu_R*l*denom*reinv;
      matrix(0, 0, 1, 1) = -mu_R/rho*reinv;
      matrix(0, 0, 2, 0) =((mu_R*ll*denom2)-hkocv*(-re*denom))*reinv;
      matrix(0, 0, 2, 1) =(-mu_R+hkocv)*l*denom*reinv;
      matrix(0, 0, 2, 2) = -hkocv/rho*reinv;
	}

	// dB/dU
	void calc_d(Tensor<double, 2> & matrix)
	{
		//		matrix = 0;
		matrix.clear();
	}
};


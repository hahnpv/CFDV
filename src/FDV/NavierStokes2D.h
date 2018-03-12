#pragma once

#include "tau.h"
#include "Flow.h"
#include "Thermo.h"
#include <cmath>
#include "NavierStokes.h" 
struct NavierStokes2D : public NavierStokes
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
		matrix(0, 2) = (l*m) / rho;
		matrix(0, 3) = (p + e)*(l/rho);

		matrix(1, 0) = m;
		matrix(1, 1) = (l*m)/rho;
		matrix(1, 2) = p + pow(m,2)/rho;
		matrix(1, 3) = (p + e)*(m/rho);
	}

	// 2D Viscous Flux
	void calc_G(Tau & tau, Tensor<double, 1> & v, Tensor<double, 1> & dT, double & HK, double & Re, Tensor<double, 2> & matrix)
	{	
		matrix.clear();

		double reinv = 1.0 / Re;

		matrix(0, 0) =0.0;
		matrix(0, 1) =-tau(0, 0)*reinv;
		matrix(0, 2) =-tau(0, 1)*reinv;
		matrix(0, 3) =(-tau(0, 0)*v(0) - tau(0, 1)*v(1) - HK*dT(0))*reinv;

		matrix(1, 0) =0.0;
		matrix(1, 1) =-tau(1,0)*reinv;
		matrix(1, 2) =-tau(1,1)*reinv;
		matrix(1, 3) =(-tau(1, 0)*v(0) - tau(1,1)*v(1) - HK*dT(1))*reinv;
	}

	/// 2D Source Term, not incorportated
	void calc_B(Node * n, Tensor<double,1> & matrix)
	{
		matrix.clear();
	};


		// 2D Flux Jacobians

	// dF/dU
	void calc_a(Flow & f, Tensor<double,3> & matrix)
	{
		Thermo & thermo = Thermo::Instance();
		matrix.clear();

      double u  = f.v(0);
      double v  = f.v(1);
      double uu = u*u;
      double uv = u*v;
      double vv = v*v;
 //     double v11p22 = uu + vv;
	  double E = f.E;

      matrix(0, 0, 1) = 1.0;
	  matrix(0, 1, 0) = thermo.gm3d2*uu+thermo.gm1d2*vv;
	  matrix(0, 1, 1) =-thermo.gamm3*u;
	  matrix(0, 1, 2) =-thermo.gamm1*v;
	  matrix(0, 1, 3) = thermo.gamm1;
      matrix(0, 2, 0) =-uv;
      matrix(0, 2, 1) = v;
      matrix(0, 2, 2) = u;
	  matrix(0, 3, 0) =-thermo.gamma*E*u+thermo.gamm1*u*(uu + vv);
	  matrix(0, 3, 1) = thermo.gamma*E-thermo.gm1d2*(3.*uu + vv);
	  matrix(0, 3, 2) =-thermo.gamm1*uv;
      matrix(0, 3, 3) = thermo.gamma*u;

	  matrix(1, 0, 2) = 1.0;
      matrix(1, 1, 0) =-uv;
      matrix(1, 1, 1) = v;
      matrix(1, 1, 2) = u;
	  matrix(1, 2, 0) = thermo.gm1d2*uu+thermo.gm3d2*vv;
	  matrix(1, 2, 1) =-thermo.gamm1*u;
	  matrix(1, 2, 2) =-thermo.gamm3*v;
	  matrix(1, 2, 3) = thermo.gamm1;
	  matrix(1, 3, 0) =-thermo.gamma*E*v+thermo.gamm1*v*(uu + vv);
	  matrix(1, 3, 1) =-thermo.gamm1*uv;
	  matrix(1, 3, 2) = thermo.gamma*E-thermo.gm1d2*(uu+3.*vv);
      matrix(1, 3, 3) = thermo.gamma*v;
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
      double v = flow.v(1);
      double uu = u*u;
      double vv = v*v;
      double v11p22 = uu + vv;

	  double drdx = flow.dU(0, 0);
      double dldx = flow.dU(1, 0);
      double dmdx = flow.dU(2, 0);
      double dedx = flow.dU(3, 0);
      double drdy = flow.dU(0, 1);
      double dldy = flow.dU(1, 1);
      double dmdy = flow.dU(2, 1);
      double dedy = flow.dU(3, 1);

		// check ordering of tau, etc.

	  matrix(0, 1, 0) =(DENOM*(mu_R*dldx+lambda*dmdy-2.*(mu_R*u*drdx+lambda*v*drdy)))*reinv;
      matrix(0, 1, 1) =mu_R*drdx*DENOM*reinv;
      matrix(0, 1, 2) =lambda*drdy*DENOM*reinv;
      matrix(0, 2, 0) =(mu*DENOM*(dldy+dmdx-2.*(u*drdy+v*drdx)))*reinv;
      matrix(0, 2, 1) =mu*drdy*DENOM*reinv;
      matrix(0, 2, 2) =mu*drdx*DENOM*reinv;

	  // may be issue with a few terms here, third term might need a minus sign?
	  matrix(0, 3, 0) =u*matrix(0, 1, 0)+v*matrix(0, 2, 0)+((u*tau(0,0)+v*tau(1,0))/rho -TERM1*(-dedx+(2.*E-3.*uu-3.*vv)*drdx+2.*u*dldx+2.*v*dmdx))*reinv;
      matrix(0, 3, 1) =u*matrix(0, 1, 1)+v*matrix(0, 2, 1)+(-tau(0,0)/rho-TERM1*(2.*u*drdx-dldx))*reinv;
      matrix(0, 3, 2) =u*matrix(0, 1, 2)+v*matrix(0, 2, 2)+(-tau(1,0)/rho-TERM1*(-dmdx+2.*v*drdx))*reinv;
      matrix(0, 3, 3) =TERM1*drdx*reinv;

	  matrix(1, 1, 0) =matrix(0, 2, 0);
      matrix(1, 1, 1) =matrix(0, 2, 1);
      matrix(1, 1, 2) =matrix(0, 2, 2);
      matrix(1, 2, 0) =(DENOM*(lambda*(dldx+4.*v*drdy)+mu_R*(u*drdx+dmdy)))*reinv;
      matrix(1, 2, 1) =lambda*drdx*DENOM*reinv;
      matrix(1, 2, 2) =mu_R*drdy*DENOM*reinv;
      matrix(1, 3, 0) =u*matrix(1, 1, 0)+v*matrix(1, 2, 0)+((u*tau(0,1)+v*tau(1,1))/rho 
              -TERM1*(-dedy+(2.*E-3.* uu-3.*vv)*drdy+2.*u*dldy+2.*v*dmdy))*reinv;
      matrix(1, 3, 1) =u*matrix(1, 1, 1)+v*matrix(1, 2, 1)+(-tau(1,0)/rho - TERM1*(-dldy+2.*u*drdy))*reinv;
      matrix(1, 3, 2) =u*matrix(1, 1, 2)+v*matrix(1, 2, 2)+(-tau(1,1)/rho - TERM1*(2.*v*drdy-dmdy))*reinv;
      matrix(1, 3, 3) =TERM1*drdy*reinv;
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
      double m     = flow.U(2);
      double re      = flow.U(3);
      double ll    = l * l;
      double lm    = l * m;
      double mm    = m * m;
      double rv11p22 = ll + mm;
	  double hkocv   = flow.k/thermo.Cv;

	  // good to here

      matrix(0, 0, 1, 0) =  mu_R*l*denom*reinv;
      matrix(0, 0, 1, 1) = -mu_R/rho*reinv;
      matrix(0, 0, 2, 0) =  mu*m*denom*reinv;
      matrix(0, 0, 2, 2) = -mu/rho*reinv;
      matrix(0, 0, 3, 0) =((mu_R*ll*denom2+mu*mm*denom2)-hkocv*(-re*denom+(rv11p22)*denom2))*reinv;
      matrix(0, 0, 3, 1) =(-mu_R+hkocv)*l*denom*reinv;
      matrix(0, 0, 3, 2) =(-mu +hkocv)*m*denom*reinv;
      matrix(0, 0, 3, 3) = -hkocv/rho*reinv;

      matrix(0, 1, 1, 0) = lambda*m*denom*reinv;
      matrix(0, 1, 1, 2) =-lambda/rho*reinv;
      matrix(0, 1, 2, 0) = mu*l*denom*reinv;
      matrix(0, 1, 2, 1) =-mu/rho*reinv;
      matrix(0, 1, 3, 0) =(lambda+mu)*lm*denom2*reinv;
      matrix(0, 1, 3, 1) =-mu*m*denom*reinv;
      matrix(0, 1, 3, 2) =-lambda*l*denom*reinv;

      matrix(1, 0, 1, 0) = mu*m*denom*reinv;
      matrix(1, 0, 1, 2) =-mu/rho*reinv;
      matrix(1, 0, 2, 0) = lambda*l*denom*reinv;
      matrix(1, 0, 2, 1) =-lambda/rho*reinv;
      matrix(1, 0, 3, 0) =(mu+lambda)*lm*denom2*reinv;
      matrix(1, 0, 3, 1) =-lambda*m*denom*reinv;
      matrix(1, 0, 3, 2) =-mu*l*denom*reinv;

      matrix(1, 1, 1, 0) = mu*l*denom*reinv;
      matrix(1, 1, 1, 1) =-mu/rho*reinv;
      matrix(1, 1, 2, 0) = mu_R*m*denom*reinv;
      matrix(1, 1, 2, 2) =-mu_R/rho*reinv;
      matrix(1, 1, 3, 0) =((mu*ll*denom2+mu_R*mm*denom2) -hkocv*(-re*denom+rv11p22*denom2))*reinv;
	  matrix(1, 1, 3, 1) =(-mu +hkocv)*l*denom*reinv;
      matrix(1, 1, 3, 2) =(-mu_R+hkocv)*m*denom*reinv;
      matrix(1, 1, 3, 3) =-hkocv/rho*reinv;
	}

	// dB/dU
	void calc_d(Tensor<double, 2> & matrix)
	{
		//		matrix = 0;
		matrix.clear();
	}
};


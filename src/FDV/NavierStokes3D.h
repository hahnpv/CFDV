#pragma once

#include "Thermo.h"
#include "NavierStokes.h"

struct NavierStokes3D : public NavierStokes
{
	/// F vector of the Navier-Stokes equation
	/// Validated 01/24/09 to Appendix A, Chung.
	void calc_F(Tensor<double, 1> & U, double & p, Tensor<double, 2> & matrix)
	{
		matrix.clear();

		double rho = U(0);
		double u = U(1)/U(0);
		double v = U(2)/U(0);
		double w = U(3)/U(0);
		double E = U(4)/U(0);	

		matrix(0, 0) = rho*u;
		matrix(0, 1) = rho*u*u + p;
		matrix(0, 2) = rho*u*v;
		matrix(0, 3) = rho*u*w;
		matrix(0, 4) = (rho*E + p)*u;

		matrix(1, 0) = rho*v;
		matrix(1, 1) = rho*u*v;
		matrix(1, 2) = rho*v*v + p;
		matrix(1, 3) = rho*v*w;
		matrix(1, 4) = (rho*E + p)*v;

		matrix(2, 0) = rho*w;
		matrix(2, 1) = rho*u*w;
		matrix(2, 2) = rho*v*w;
		matrix(2, 3) = rho*w*w + p;
		matrix(2, 4) = (rho*E + p)*w;
	}

	/// G vector of the Navier-Stokes equation
	/// Validated 01/24/09 to Appendix A, Chung.
	void calc_G(Tau & tau, Tensor<double, 1> & v, Tensor<double, 1> & dT, double & HK, double & Re, Tensor<double, 2> & matrix)
	{
		matrix.clear();

		double reinv = 1.0 / Re;

		matrix(0, 0) = 0.0;
		matrix(0, 1) = -tau(0,0)*reinv;
		matrix(0, 2) = -tau(0,1)*reinv;
		matrix(0, 3) = -tau(0,2)*reinv;
		matrix(0, 4) = (-tau(0,0)*v(0) - tau(0,1)*v(1) -tau(0,2)*v(2) - HK*dT(0))*reinv;

		matrix(1, 0) = 0.0;
		matrix(1, 1) = -tau(1,0)*reinv;
		matrix(1, 2) = -tau(1,1)*reinv;
		matrix(1, 3) = -tau(1,2)*reinv;
		matrix(1, 4) = (-tau(1,0)*v(0) - tau(1,1)*v(1) -tau(1,2)*v(2) - HK*dT(1))*reinv;

		matrix(2, 0) = 0.0;
		matrix(2, 1) = -tau(2,0)*reinv;
		matrix(2, 2) = -tau(2,1)*reinv;
		matrix(2, 3) = -tau(2,2)*reinv;
		matrix(2, 4) = (-tau(2,0)*v(0) - tau(2,1)*v(1) -tau(2,2)*v(2) - HK*dT(2))*reinv;
	}

	// Source Term
	void calc_B(Node * n, Tensor<double,1> & matrix)
	{
		matrix.clear();
	};

		/// dF/dU
		/// Validated 01/24/09 to Appendix A, Chung.
	void calc_a(Flow & f, Tensor<double,3> & matrix)
	{
		Thermo & thermo = Thermo::Instance();
		double gamma = thermo.gamma;
		double gamma1 =	thermo.gamm1;
		double gamma3 =	thermo.gamm3;
		double gm3d2 = thermo.gm3d2;
		double gamm3 = thermo.gamm3;
		double gm1d2 = thermo.gm1d2;

		matrix.clear();

      double u  = f.v(0);
      double v  = f.v(1);
	  double w  = f.v(2);
      double uu = u*u;
      double uv = u*v;
      double vv = v*v;
	  double ww = w*w;
	  double uw = u*w;
	  double vw = v*w;
	  double E = f.E;		// U(4) / U(0)

	  matrix(0, 0, 0) = 0.0;
      matrix(0, 0, 1) = 1.0;
	  matrix(0, 0, 2) = 0.0;
	  matrix(0, 0, 3) = 0.0;
	  matrix(0, 0, 4) = 0.0;

      matrix(0, 1, 0)=gm3d2*uu+gm1d2*(vv + ww);
      matrix(0, 1, 1)=-gamma3*u;
      matrix(0, 1, 2)=-gamma1*v;
      matrix(0, 1, 3)=-gamma1*w;
      matrix(0, 1, 4)=gamma1;

      matrix(0, 2, 0)=-uv;
      matrix(0, 2, 1)=v;
      matrix(0, 2, 2)=u;
	  matrix(0, 2, 3) = 0.0;
	  matrix(0, 2, 4) = 0.0;

      matrix(0, 3, 0) =-uw;
      matrix(0, 3, 1) = w;
	  matrix(0, 3, 2) = 0.0;
      matrix(0, 3, 3) = u;
	  matrix(0, 3, 4) = 0.0;

      matrix(0, 4, 0)=-gamma*E*u+gamma1*u*(uu + vv + ww);
      matrix(0, 4, 1)=gamma*E-gm1d2*(3.*uu+(vv + ww));
      matrix(0, 4, 2)=-gamma1*uv;
      matrix(0, 4, 3)=-gamma1*uw;
      matrix(0, 4, 4)=gamma*u;


	  matrix(1, 0, 0) = 0.0;
	  matrix(1, 0, 1) = 0.0;
      matrix(1, 0, 2) = 1.0;
	  matrix(1, 0, 3) = 0.0;
	  matrix(1, 0, 4) = 0.0;

      matrix(1, 1, 0) =-uv;
      matrix(1, 1, 1) = v;
      matrix(1, 1, 2) = u;
	  matrix(1, 1, 3) = 0.0;
	  matrix(1, 1, 4) = 0.0;

      matrix(1, 2, 0)=gm1d2*(uu + ww)+gm3d2*vv;
      matrix(1, 2, 1)=-gamma1*u;
      matrix(1, 2, 2)=-gamma3*v;
      matrix(1, 2, 3)=-gamma1*w;
      matrix(1, 2, 4)=gamma1;

      matrix(1, 3, 0) = -vw;
 	  matrix(1, 3, 1) = 0.0;
      matrix(1, 3, 2) = w;
      matrix(1, 3, 3) = v;
	  matrix(1, 3, 4) = 0.0;

      matrix(1, 4, 0)=-gamma*E*v+gamma1*v*(uu + vv + ww);
      matrix(1, 4, 1)=-gamma1*uv;
      matrix(1, 4, 2)=gamma*E-gm1d2*((uu + ww)+3.*vv);
      matrix(1, 4, 3)=-gamma1*vw;
      matrix(1, 4, 4)=gamma*v;


      matrix(2, 0, 0) = 0.0;
	  matrix(2, 0, 1) = 0.0;
	  matrix(2, 0, 2) = 0.0;
	  matrix(2, 0, 3) = 1.0;
  	  matrix(2, 0, 4) = 0.0;

      matrix(2, 1, 0) =-uw;
      matrix(2, 1, 1) = w;
  	  matrix(2, 1, 2) = 0.0;
      matrix(2, 1, 3) = u;
	  matrix(2, 1, 4) = 0.0;


      matrix(2, 2, 0) = -vw;
  	  matrix(2, 2, 1) = 0.0;
      matrix(2, 2, 2) = w;
      matrix(2, 2, 3) = v;
  	  matrix(2, 2, 4) = 0.0;

      matrix(2, 3, 0)=gm1d2*(uu + vv) + gm3d2*ww;
      matrix(2, 3, 1)=-gamma1*u;
      matrix(2, 3, 2)=-gamma1*v;
      matrix(2, 3, 3)=-gamma3*w;
      matrix(2, 3, 4)=gamma1;

      matrix(2, 4, 0)=-gamma*E*w + gamma1*w*(uu + vv + ww);
      matrix(2, 4, 1)=-gamma1*uw;
      matrix(2, 4, 2)=-gamma1*vw;
      matrix(2, 4, 3)=gamma*E - gm1d2*((uu + vv) + 3.*ww);
      matrix(2, 4, 4)=gamma*w;
	}

		/// dG/dU
		/// Validated 01/24/09 to Appendix A, Chung.
	void calc_b(Flow & flow, Tensor<double, 3> & matrix)
	{
		Thermo & thermo = Thermo::Instance();

		matrix.clear();

//		double reinv = 1.0/thermo.creyn;
		double reinv = 1.0;
		double mu = flow.mu;

		Tensor<double, 2> &tau = flow.tau;			// try without reference now that copyop and assign should be fixed.

		double E = flow.E;

		double lambda=-2.*mu/3.;
		double mu_R=2.*mu+lambda;
		double denom=1./pow(flow.U(0),2);
		double TERM1=flow.k*denom/thermo.Cv;

		double rho = flow.U(0);
		double u = flow.v(0);
		double v = flow.v(1);
		double w = flow.v(2);
		double uu = u*u;
		double vv = v*v;
		double ww = w*w;
		double v11p22 = uu + vv;

		double drdx = flow.dU(0, 0);		// verify order
		double dldx = flow.dU(1, 0);
		double dmdx = flow.dU(2, 0);
		double dndx = flow.dU(3, 0);
		double dedx = flow.dU(4, 0);

		double drdy = flow.dU(0, 1);
		double dldy = flow.dU(1, 1);
		double dmdy = flow.dU(2, 1);
		double dndy = flow.dU(3, 1);
		double dedy = flow.dU(4, 1);

		double drdz = flow.dU(0, 2);
		double dldz = flow.dU(1, 2);
		double dmdz = flow.dU(2, 2);
		double dndz = flow.dU(3, 2);
		double dedz = flow.dU(4, 2);

		double l = flow.U(1);		// ru
		double m = flow.U(2);		// rv
		double n = flow.U(3);		// rw

      matrix(0, 1, 0)=(denom*(mu_R*dldx+lambda*(dmdy+dndz)+ mu_R*(-2.*u*drdx+v*drdy+ w*drdz)))*reinv;
      matrix(0, 1, 1)=mu_R*drdx*denom*reinv;
      matrix(0, 1, 2)=lambda*drdy*denom*reinv;
      matrix(0, 1, 3)=lambda*drdz*denom*reinv;
      matrix(0, 2, 0)=(mu*denom*(dldy+dmdx-2.*(u*drdy+v*drdx)))*reinv;
      matrix(0, 2, 1)=mu*drdy*denom*reinv;
      matrix(0, 2, 2)=mu*drdx*denom*reinv;
      matrix(0, 3, 0)=(mu*denom*(dldz+dndx-2.*(u*drdz+w*drdx)))*reinv;
      matrix(0, 3, 1)=mu*drdz*denom*reinv;
      matrix(0, 3, 3)=mu*drdx*denom*reinv;
      matrix(0, 4, 0)=u*matrix(0, 1, 0)+v*matrix(0, 2, 0)+w*matrix(0, 3, 0)+(denom*(l*tau(0,0) +
              m*tau(0,1)+n*tau(0,2))
			  -TERM1*(-dedx+(2.*E-3.*(uu + vv + ww))*drdx    
			  +2.*(u*dldx+v*dmdx+w*dndx)))*reinv;

	  // done to here.
/*
      matrix(0, 4, 0)=u*matrix(0, 1, 0)+v*matrix(0, 2, 0)+w*matrix(0, 3, 0)+(denom*(l*tau(0,0) +
              m*tau(0,1)+n*tau(0,2))			// changed order of tau indices per gwh dissertation text pg 144
			  -TERM1*(-dedx+(2.*E-3.*(uu + vv + ww))*drdx    
			  +2.*(u*dldx+v*dmdx+w*dndx)))*reinv;
*/
	  matrix(0, 4, 1)=u*matrix(0, 1, 1)+v*matrix(0, 2, 1)+w*matrix(0, 3, 1) + (
              -tau(0,0)/rho - TERM1*(2.*u*drdx-dldx))*reinv;

      matrix(0, 4, 2)=u*matrix(0, 1, 2)+v*matrix(0, 2, 2)+w*matrix(0, 3, 2) + (-tau(0,1)/rho   
              -TERM1*(-dmdx+2.*v*drdx))*reinv;

      matrix(0, 4, 3)=u*matrix(0, 1, 3)+v*matrix(0, 2, 3)+w*matrix(0, 3, 3) + (-tau(0,2)/rho - 
              TERM1*(2*w*drdx-dndx))*reinv;
      matrix(0, 4, 4)=TERM1*drdx*reinv;

      matrix(1, 1, 0)=matrix(0, 2, 0);
      matrix(1, 1, 1)=matrix(0, 2, 1);
      matrix(1, 1, 2)=matrix(0, 2, 2);

      matrix(1, 2, 0)=(denom*(lambda*(dldx+dndz) + mu_R*(dmdy+u*drdx-2.*v*  
              drdy+w*drdz)))*reinv;
	  matrix(1, 2, 1)=lambda*drdx*denom*reinv;
      matrix(1, 2, 2)=mu_R*drdy*denom*reinv;
      matrix(1, 2, 3)=lambda*drdz*denom*reinv;
      matrix(1, 3, 0)=(mu*denom*(dmdz+dndy-2.*(v*drdz +w*drdy)))*reinv;
      matrix(1, 3, 2)=mu*denom*drdz*reinv;
      matrix(1, 3, 3)=mu*denom*drdy*reinv;

      matrix(1, 4, 0)=u*matrix(1, 1, 0)+v*matrix(1, 2, 0)+w*matrix(1, 3, 0) + (denom*(l*
              tau(1,0)+m*tau(1,1)+n*tau(1,2))		// change per dissertation (one typo here)
              -TERM1*(-dedy+(2.*E-3.*(uu + vv + ww))	
              *drdy + 2.*(u*dldy+v*dmdy+w*dndy)))*reinv;

/*
matrix(1, 4, 0)=u*matrix(1, 1, 0)+v*matrix(1, 2, 0)+w*matrix(1, 3, 0) + (denom*(l*
              tau(1,0)+m*tau(1,1)+n*tau(1,2))		// change per dissertation (one typo here)
              -TERM1*(-dedy+(2.*E-3.*(uu + vv + ww))	
              *drdy + 2.*(u*dldy+v*dmdy+w*dndy)))*reinv;
*/
      matrix(1, 4, 1)=u*matrix(1, 1, 1) + v*matrix(1, 2, 1) + w*matrix(1, 3, 1)+(-tau(1,0)/rho - 
              TERM1*(-dldy+2.*u*drdy))*reinv;
      matrix(1, 4, 2)=u*matrix(1, 1, 2)+v*matrix(1, 2, 2)+w*matrix(1, 3, 2)+            
              (-tau(1,1)/rho - TERM1*(2.*v*drdy - dmdy))*reinv;
      matrix(1, 4, 3)=u*matrix(1, 1, 3)+v*matrix(1, 2, 3)+w*matrix(1, 3, 3) + (-tau(1,2)/rho    
              -TERM1*(2.*w*drdy-dndy))*reinv;
      matrix(1, 4, 4)=TERM1*drdy*reinv;

      matrix(2, 1, 0)=matrix(0, 3, 0);
      matrix(2, 1, 1)=matrix(0, 3, 1);
      matrix(2, 1, 3)=matrix(0, 3, 3);
      matrix(2, 2, 0)=matrix(1, 3, 0);
      matrix(2, 2, 2)=matrix(1, 3, 2);
      matrix(2, 2, 3)=matrix(1, 3, 3);

      matrix(2, 3, 0)=(denom*(mu_R*(dndz+u*drdx+v*drdy-2.*w*drdz)+lambda*(dldx+ dmdy)))*reinv;
      matrix(2, 3, 1)=lambda*drdx*denom*reinv;
      matrix(2, 3, 2)=lambda*drdy*denom*reinv;
      matrix(2, 3, 3)=mu_R*drdz*denom*reinv;
      matrix(2, 4, 0)=u*matrix(2, 1, 0)+v*matrix(2, 2, 0)+w*matrix(2, 3, 0) + (
              denom*(l*tau(2,0)+m*tau(2,1)+n*	
              tau(2,2))-TERM1*(-dedz+(2.*E-3.*(uu + vv + ww))*drdz  
              +2.*(u*dldz+v*dmdz+w*dndz)))*reinv;
      matrix(2, 4, 1)=u*matrix(2, 1, 1)+v*matrix(2, 2, 1)+w*matrix(2, 3, 1) + (-tau(2,0)/rho  
              -TERM1*(2.*u*drdz-dldz))*reinv;
      matrix(2, 4, 2)=u*matrix(2, 1, 2)+v*matrix(2, 2, 2)+w*matrix(2, 3, 2) + (-tau(2,1)/rho  
              -TERM1*(2.*v*drdz-dmdz))*reinv;
      matrix(2, 4, 3)=u*matrix(2, 1, 3)+v*matrix(2, 2, 3)+w*matrix(2, 3, 3) +               
              (-tau(2,2)/rho - TERM1*(2.*w*drdz - dndz))*reinv;
      matrix(2, 4, 4)=TERM1*drdz*reinv;


	  // test
	  reinv = 1.0/thermo.creyn;
	  for (int i = 0; i < matrix.imax(); i++)
	  {
		  for (int j = 0; j < matrix.jmax(); j++)
		  {
			  for (int k = 0; k < matrix.kmax(); k++)
			  {
					matrix(i, j, k) *= reinv;
			  }
		  }
	  }
	}

	// dG/dU,j
	/// Validated 01/24/09 to Appendix A, Chung.
	void calc_c(Flow & flow, Tensor<double, 4> & matrix)
	{
		Thermo & thermo = Thermo::Instance();

		matrix.clear();

//		double reinv = 1.0/thermo.creyn;
		double reinv = 1.0;
		double mu = flow.mu;

      double lambda=-2.*mu/3.;
      double mu_R=2.*mu+lambda;

      double rho      = flow.U(0);
      double denom  = 1./pow(rho,2);

      double l     = flow.U(1);			// ru
      double m     = flow.U(2);			// rv
	  double n		= flow.U(3);		// rw
      double re      = flow.U(4);		// rho * E = epsilon

      double ll    = l * l;
	  double mm    = m * m;
	  double nn		= n * n;
	  double hkocv   = flow.k/thermo.Cv;

	  double TERM1 = flow.k/thermo.Cv * denom;

      matrix(0, 0, 1, 0)=mu_R*l*denom*reinv;
      matrix(0, 0, 1, 1)=(-mu_R/rho)*reinv;
      matrix(0, 0, 2, 0)=mu*m*denom*reinv;
      matrix(0, 0, 2, 2)=(-mu/rho)*reinv;
      matrix(0, 0, 3, 0)=mu*n*denom*reinv;
      matrix(0, 0, 3, 3)=(-mu/rho)*reinv;
      matrix(0, 0, 4, 0)=(l*matrix(0, 0, 1, 0)+ m*matrix(0, 0, 2, 0) + n*matrix(0, 0, 3, 0))/rho
                -TERM1*(-re + (ll + mm + nn)/rho)*reinv;
      matrix(0, 0, 4, 1)=(l/rho)*matrix(0, 0, 1, 1) + TERM1*l*reinv;
      matrix(0, 0, 4, 2)=(m/rho)*matrix(0, 0, 2, 2) + TERM1*m*reinv;
      matrix(0, 0, 4, 3)=(n/rho)*matrix(0, 0, 3, 3) + TERM1*n*reinv;
      matrix(0, 0, 4, 4)=(-hkocv/rho)*reinv;

      matrix(0, 1, 1, 0)=lambda*m*denom*reinv;
      matrix(0, 1, 1, 2)=(-lambda/rho)*reinv;
      matrix(0, 1, 2, 0)=mu*l*denom*reinv;
      matrix(0, 1, 2, 1)=(-mu/rho)*reinv;
      matrix(0, 1, 4, 0)=(l/rho)*matrix(0, 1, 1, 0)+(m/rho)*matrix(0, 1, 2, 0);
      matrix(0, 1, 4, 1)=(m/rho)*matrix(0, 1, 2, 1);
      matrix(0, 1, 4, 2)=(l/rho)*matrix(0, 1, 1, 2);

      matrix(0, 2, 1, 0)=lambda*n*denom*reinv;
      matrix(0, 2, 1, 3)=(-lambda/rho)*reinv;
      matrix(0, 2, 3, 0)=mu*l*denom*reinv;
      matrix(0, 2, 3, 1)=(-mu/rho)*reinv;
      matrix(0, 2, 4, 0)=(l/rho)*matrix(0, 2, 1, 0)+(n/rho)*matrix(0, 2, 3, 0);
      matrix(0, 2, 4, 1)=(n/rho)*matrix(0, 2, 3, 1);
      matrix(0, 2, 4, 3)=(l/rho)*matrix(0, 2, 1, 3);

      matrix(1, 0, 1, 0)=matrix(0, 0, 2, 0);
      matrix(1, 0, 1, 2)=matrix(0, 0, 2, 2);
      matrix(1, 0, 2, 0)=lambda*l*denom*reinv;
      matrix(1, 0, 2, 1)=(-lambda/rho)*reinv;
      matrix(1, 0, 4, 0)=(l/rho)*matrix(1, 0, 1, 0)+(m/rho)*matrix(1, 0, 2, 0);
      matrix(1, 0, 4, 1)=(m/rho)*matrix(1, 0, 2, 1);
      matrix(1, 0, 4, 2)=(l/rho)*matrix(1, 0, 1, 2);

      matrix(1, 1, 1, 0)=matrix(0, 1, 2, 0);
      matrix(1, 1, 1, 1)=matrix(0, 1, 2, 1);
      matrix(1, 1, 2, 0)=mu_R*m*denom*reinv;
      matrix(1, 1, 2, 2)=(-mu_R/rho)*reinv;
      matrix(1, 1, 3, 0)=mu*n*denom*reinv;
      matrix(1, 1, 3, 3)=(-mu/rho)*reinv;
      matrix(1, 1, 4, 0)=(l*matrix(1, 1, 1, 0)+m*matrix(1, 1, 2, 0)+n*matrix(1, 1, 3, 0))/rho 
                - TERM1*(-re + (ll + mm + nn)/rho)*reinv;
      matrix(1, 1, 4, 1)=(l/rho)*matrix(1, 1, 1, 1) + TERM1*l*reinv;
      matrix(1, 1, 4, 2)=(m/rho)*matrix(1, 1, 2, 2) + TERM1*m*reinv;
      matrix(1, 1, 4, 3)=(n/rho)*matrix(1, 1, 3, 3) + TERM1*n*reinv;
      matrix(1, 1, 4, 4)=(-hkocv/rho)*reinv;

      matrix(1, 2, 2, 0)=lambda*n*denom*reinv;
      matrix(1, 2, 2, 3)=(-lambda/rho)*reinv;
      matrix(1, 2, 3, 0)=mu*m*denom*reinv;
      matrix(1, 2, 3, 2)=(-mu/rho)*reinv;
      matrix(1, 2, 4, 0)=(m/rho)*matrix(1, 2, 2, 0)+(n/rho)*matrix(1, 2, 3, 0);
      matrix(1, 2, 4, 2)=(n/rho)*matrix(1, 2, 3, 2);
      matrix(1, 2, 4, 3)=(m/rho)*matrix(1, 2, 2, 3);

      matrix(2, 0, 1, 0)=matrix(0, 0, 3, 0);
      matrix(2, 0, 1, 3)=matrix(0, 0, 3, 3);
      matrix(2, 0, 3, 0)=lambda*l*denom*reinv;
      matrix(2, 0, 3, 1)=(-lambda/rho)*reinv;
      matrix(2, 0, 4, 0)=(l/rho)*matrix(2, 0, 1, 0)+(n/rho)*matrix(2, 0, 3, 0);
      matrix(2, 0, 4, 1)=(n/rho)*matrix(2, 0, 3, 1);
      matrix(2, 0, 4, 3)=(l/rho)*matrix(2, 0, 1, 3);
      matrix(2, 1, 2, 0)=matrix(1, 1, 3, 0);
      matrix(2, 1, 2, 3)=matrix(1, 1, 3, 3);
      matrix(2, 1, 3, 0)=lambda*m*denom*reinv;
      matrix(2, 1, 3, 2)=(-lambda/rho)*reinv;
      matrix(2, 1, 4, 0)=(m/rho)*matrix(2, 1, 2, 0)+(n/rho)*matrix(2, 1, 3, 0);
      matrix(2, 1, 4, 2)=(n/rho)*matrix(2, 1, 3, 2);
      matrix(2, 1, 4, 3)=(m/rho)*matrix(2, 1, 2, 3);

      matrix(2, 2, 1, 0)=matrix(0, 2, 3, 0);
      matrix(2, 2, 1, 1)=matrix(0, 2, 3, 1);
      matrix(2, 2, 2, 0)=matrix(1, 2, 3, 0);
      matrix(2, 2, 2, 2)=matrix(1, 2, 3, 2);
      matrix(2, 2, 3, 0)=mu_R*n*denom*reinv;
      matrix(2, 2, 3, 3)=(-mu_R/rho)*reinv;
      matrix(2, 2, 4, 0)=(l*matrix(2, 2, 1, 0)+m*matrix(2, 2, 2, 0)+n*matrix(2, 2, 3, 0))/rho
                - TERM1*(-re + (ll + mm + nn)/rho)*reinv;
      matrix(2, 2, 4, 1)=(l/rho)*matrix(2, 2, 1, 1) + TERM1*l*reinv;
      matrix(2, 2, 4, 2)=(m/rho)*matrix(2, 2, 2, 2) + TERM1*m*reinv;
      matrix(2, 2, 4, 3)=(n/rho)*matrix(2, 2, 3, 3) + TERM1*n*reinv;
      matrix(2, 2, 4, 4)=(-hkocv/rho)*reinv;


	  	  // test
	  reinv = 1.0/thermo.creyn;
	  for (int i = 0; i < matrix.imax(); i++)
	  {
		  for (int j = 0; j < matrix.jmax(); j++)
		  {
			  for (int k = 0; k < matrix.kmax(); k++)
			  {
				  for (int l = 0; l < matrix.lmax(); l++)
				  {
					matrix(i, j, k, l) *= reinv;
				  }
			  }
		  }
	  }
	}

	// dB/dU
	void calc_d(Tensor<double, 2> & matrix)
	{
		matrix.clear();
	}

};


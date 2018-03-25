#pragma once

#include <iostream>
#include <numeric>
using namespace std;

	// 1D
#include "NavierStokes1D.h"
#include "../FiniteElement/TestFunctions1D.h"

	// 2D
#include "NavierStokes2D.h"
#include "Function/eval1dboundary.h"
#include "../FiniteElement/TestFunctions1D.h"

	// 3D
#include "NavierStokes3D.h"
#include "../FiniteElement/calc_TF_norm3D.h"					// norm + test functions
#include "Function/eval2dboundary.h"

	// Quadrature
#include "../FiniteElement/GaussQuad.h"

#include "Function/FDVParam.h"
	// Adaptive Meshing
#include "../adap/MeshRefine2D.h"

#include "Flow.h"

// temporary
// #include "PrecisionTimer.h"

	// This is the Finite Element implementation of the FDV method as presented by Dr. T. J. Chung
	// in "Computational Fluid Dynamics", chapter 13, using the standard Galerkin method.

	// TODO: 
	// 1. Jacobian swaps 2d/3d
	// 2. F/G swaps in node 2d/3d

	/// Solve the FDV system of equations formulated using the standard Galerkin method.
// : g and k, h and m, p and t, q and u. are the opposite sign of each other, former used in domain,
// latter in the boundary. If we flip signs in the boundary formation, instead of keeping the signage in the matrices,
// we could probably get away with only
// using g, h, p, q summations and have a "unified" matrix summation.

double coth( double x)
{
	return cosh(x) / sinh(x);

}

template<class T >
struct FDVGalerkin : public unary_function<T, void>
{
	FDVGalerkin(const int nnod, const int neqn, const int ndim, const int nbnod, double dt,NavierStokes * NS)
		:	g(ndim,neqn,neqn),					
			h(ndim,ndim,neqn,neqn),
			p(nnod,ndim,neqn),
			q(nnod,ndim,ndim,neqn),
			eta(neqn,neqn),
			t(nbnod,ndim,neqn),
			u(nbnod,ndim,neqn,ndim),
			k(ndim,neqn,neqn),
			m(ndim,ndim,neqn,neqn),
			s(nnod,neqn),
			flow(ndim,neqn),
			ndim(ndim),
			neqn(neqn),
			nnod(nnod),
			nbnod(nbnod),				// number of boundary nodes
			gq(2),	
			dt(dt),
			dtsq_2(pow(dt,2.0)/2.0),
			a(ndim,neqn,neqn),
			b(ndim,neqn,neqn),
			c(ndim,ndim,neqn,neqn),
			d(neqn,neqn),
			nd(ndim, neqn, neqn),
		    NS(NS)
	{
		B.resize(nnod);
		G.resize(nnod);
		F.resize(nnod);

//		NS = new NavierStokes();
		for (int i=0; i<nnod; i++)
		{
			F[i] = new Tensor<double, 2>(ndim, neqn);
			G[i] = new Tensor<double, 2>(ndim, neqn);
			B[i] = new Tensor<double, 1>(neqn);
		}
	}
/*
	FDVGalerkin(const FDVGalerkin<T, NavierStokes> & rhs)									// CopyOp
	{
		cout << "CopyOp" << endl;
//		ref = rhs.ref;
		NS = rhs.NS;
	}
*/
	~FDVGalerkin()
	{
//		delete NS;

		for (int i=0; i<nnod; i++)
		{
			delete F[i];
			delete G[i];
			delete B[i];
		}		
	}

	FDVGalerkin<T > & operator[](double deltat)
	{
		dt = deltat;
		dtsq_2 = pow(deltat, 2.0) / 2.0;
		return *this;
	}

	typedef void result_type;					// result_type for tr1 wrapper.
	result_type operator() (T e) 
	{
		element = e;

		FDVParam<Element *> fdvparam;			/// Calculate the FDV parameters

		fdvparam(e, s1, s2, s3, s4);
		s5 = 0.0;
		s6 = 0.0;

		// 2D Element Quadrature
		if (ndim == 2)
		{
			gq.two(this,&FDVGalerkin<T>::domain);
			gq.one(this,&FDVGalerkin<T>::boundary);
		}

		// 3D Element Quadrature
		else if (ndim == 3)
		{
			gq.three(this,&FDVGalerkin<T>::domain);
			gq.two(this,&FDVGalerkin<T>::boundary);
		}
		else if (ndim == 1)
		{
			gq.one(this,&FDVGalerkin<T>::domain);
//			gq.one(this,&FDVGalerkin<T, NavierStokes>::boundary);			// not sure how to handle ... 
			Tensor<double, 1> w(1);
			w(0) = 1.;
			Tensor<double, 1> zeta(1);
			zeta(0) = 1.;
			boundary(w, zeta);
		}
	}

		/// Calculate the domain integral
	void domain(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
//		cout << "calculating domain integral" << endl;
		std::vector<Node *> & node = element->node;

		double J;
		Tensor<double, 1> Phi(nnod);
		Tensor<double, 2> dPhi(nnod, ndim);

		if (ndim == 2 && nnod == 3)	// triangle
		{
			Phi = testfunctiontri(zeta);

			dtestfunctiontri dPxi(zeta);
	
			jacobiantri jac(dPxi, node);
			dPhi = dtestfunctiondxtri(dPxi,jac);

			J = jac.det();
		}
		else if (ndim == 2 && nnod == 4)
		{
			Phi = testfunction(zeta);

			dtestfunction dPxi(zeta);

			if (element->adap)
			{
				adap(Phi, dPxi, element->Hnm);
			}
	
			jacobian jac(dPxi, node);
			dPhi = dtestfunctiondx(dPxi,jac);

			J = jac.det();
		}
		else if (ndim == 3)		// was using tf3d but need to pass Hnm, adap flag for adaption, gets messy, think it through.
		{	
			Phi = testfunction3d(zeta);
			dtestfunction3d dPxi(zeta);

			if (element->adap)
			{
				adap(Phi, dPxi, element->Hnm);
			}

			jacobian3d jac(dPxi, node);
			dPhi = dtestfunctiondx3d(dPxi,jac);

			J = jac.det();
		}
		else if (ndim == 1)
		{
			Phi = Phibou(zeta);
//			dPxi dpxi;
//			en = Normal(w,dpxi,node);
			dPhi = dPidx(*node[0],*node[1]);
			J = dPhi(0,0) * node[0]->x(0) + dPhi(1,0) * node[1]->x(0);	// not sure.
		}
		
		// 2. Calculate flow variables, viscocity, k, tau
		sum(node,Phi,dPhi);
		flow.tau.calculate(flow.mu,flow.dV);

			// NUMERICAL DIFFUSION TERM // FIXME
		Tensor<double,1> Psi(nnod);		/// numerical diffusion test function
		double tau = 0;
		
		if (ndim == 2)					// eventually generalize then have a n.d. flag.
		{
			double J_xi, J_eta;
			Tensor<double, 1> e_xi(ndim), e_eta(ndim);
			double v_xi, v_eta;
			double h_xi = 0.0;
			double h_eta = 0.0;
			double R_xi, R_eta;
			double alpha_xi, alpha_eta;
			double S;

//			jacobian jac(dtestfunction(zeta), node);
			Tensor<double, 2> dtf = dtestfunction(zeta);
			if (nnod == 4)										/// fixme using inheritance or something
			{
				jacobian jac(dtf, node);
				J_xi = pow(jac(0, 0), 2) + pow(jac(1, 0), 2);
				J_eta = pow(jac(0, 1), 2) + pow(jac(1, 1), 2);

				for (int i = 0; i < ndim; i++)
				{
					e_xi(i) = 1. / pow(J_xi, 0.5)  * jac(i, 0);	// worked as dPxi()
					e_eta(i) = 1. / pow(J_eta, 0.5) * jac(i, 1);
				}
			}
			else if (nnod == 3)
			{
				jacobiantri jac(dtf, node);
				J_xi = pow(jac(0, 0), 2) + pow(jac(1, 0), 2);
				J_eta = pow(jac(0, 1), 2) + pow(jac(1, 1), 2);

				for (int i = 0; i < ndim; i++)
				{
					e_xi(i) = 1. / pow(J_xi, 0.5)  * jac(i, 0);	// worked as dPxi()
					e_eta(i) = 1. / pow(J_eta, 0.5) * jac(i, 1);
				}
			}
			/*
			else if (nnod = 6)
			{
				jacobian3d jac(dtf, node);
				J_xi = pow(jac(0, 0), 2) + pow(jac(1, 0), 2);
				J_eta = pow(jac(0, 1), 2) + pow(jac(1, 1), 2);

				for (int i = 0; i < ndim; i++)
				{
					e_xi(i) = 1. / pow(J_xi, 0.5)  * jac(i, 0);	// worked as dPxi()
					e_eta(i) = 1. / pow(J_eta, 0.5) * jac(i, 1);
				}
			}
			*/

			v_xi  = dot(flow.v, e_xi);
			v_eta = dot(flow.v, e_eta);

			S = dot(flow.v, flow.v);

			// can generalize using test function to interpolate...  FIXME should be testfunction, testfunctiontri, ... 
			// FIXME 2d specific 
			Tensor<double, 1> xloc(2);
			xloc(0) = -1.; xloc(1) =  0.;		// left
			Tensor<double, 1> xleft = testfunction(xloc);
			xloc(0) =  1.; xloc(1) =  0.;		// right
			Tensor<double, 1> xright = testfunction(xloc);
			xloc(0) =  0.; xloc(1) = -1.;		// bottom
			Tensor<double, 1> xtop = testfunction(xloc);
			xloc(0) =  0.; xloc(1) =  1.;		// top
			Tensor<double, 1> xbottom = testfunction(xloc);
			for( int N = 0; N < nnod; N++)
			{
				h_xi  += xright(N) * node[N]->x(0) - xleft(N) * node[N]->x(0);
				h_eta += xtop(N) * node[N]->x(1)   - xbottom(N) * node[N]->x(1); 
			}

			R_xi  = flow.U(0) * v_xi  * h_xi  / flow.mu;
			R_eta = flow.U(0) * v_eta * h_eta / flow.mu;

		//	cout << "R: " << R_xi << " " << R_eta << endl;
		//	cout << "S: " << S << endl;

			// If R == 0, then alpha is invalid (division by zero)
			// Set alpha = 0
			// This is justifiable since in tau, alpha is multiplied by velocity (which is zero) and nulls it anyways.
			if (R_xi == 0)
				alpha_xi = 0;
			else
				alpha_xi = coth(R_xi / 2.) - 2. / R_xi;

			if (R_eta == 0)
				alpha_eta = 0;
			else
				alpha_eta = coth(R_eta / 2.) - 2./R_eta;
	
			if (S == 0)			// FIXME hack
				tau = 0;
			else
				tau = (1./4.) * (alpha_xi * h_xi * v_xi + alpha_eta * h_eta * v_eta) / S;

			for (int N = 0; N < nnod; N++)
			{
				for (int i = 0; i < ndim; i++)
				{
					Psi(N) += tau * flow.v(i) * dPhi(N, i);
				}
			}
		//	cout << "Psi: " << Psi << endl;

		}
			// NUMERICAL DIFFUSION TERM //
		
			// DISCONTINUITY CAPTURING // FIXME this was commented out... verify correctness
		Tensor<double, 1> A(ndim);
		Tensor<double, 1> v(ndim);			// v_i^b
		double gamma = 0;
		double tau_b = 0;

		for (int k = 0; k < ndim; k++)
		{
			for (int j = 0; j < ndim; j++)
			{
				A(j) += flow.dV(j, k) * flow.v(k);
			}
		}

		gamma = mag(flow.v) * mag(A);

		for( int i = 0; i < ndim; i++)
		{
			for ( int j = 0; j < ndim; j++)
			{
				v(i) = flow.v(i) * flow.v(j) * A(j) / gamma;
			}
		}

//		double dbg = 0;
		for (int j = 0; j < ndim; j++)
		{
			tau_b += tau * flow.v(j) * A(j) / gamma;
		//	dbg += flow.v(j) * A(j) / gamma;
		}

		if ( (tau_b - tau) < 0)
		{
//			cout << "element:    " << element->number << endl;
//			cout << "taub < tau: " << tau << " " << tau_b <<  endl;
//			cout << "AjVj/Gamma: " << dbg << " A: " << A << endl;
//			cout << "v:          " << flow.v << endl;
//			cout << "vj:         " << flow.dV << endl;
			tau_b = tau;
		}
		if ( gamma == 0)
		{
			// either v = 0 or dV = 0 -> A = 0, gamma = 0; division by 0 -1#IND
			tau_b = 0;
		}

		for (int N = 0; N < nnod; N++)
		{
			for (int i = 0; i < ndim; i++)
			{
				Psi(N) += tau_b * flow.v(i) * dPhi(N, i);
			}
		}
//		cout << "Psi: " << Psi << " tau_b: " << tau_b << " gamma: " << gamma << endl;
			// DISCONTINUITY CAPTURING //

		// 4. Calculate a, b, c jacobians and F,G fluxes
		NS->calc_a(flow, a);
		NS->calc_b(flow, b);
		NS->calc_c(flow, c);
		NS->calc_d(d);

		sumFGB(node, Phi, dPhi);

		// 5. Assemble vectors 
		sumDomainTerms();

		sumDomainMatrix(w, Phi, dPhi, Psi, J);
	}

	void sumDomainMatrix(Tensor<double, 1> & w, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi, Tensor<double, 1> & Psi, double & J)
	{
		double w_i = aggregate(w);

			// RHS
		for (int N = 0; N < nnod; N++)
		{
			for (int ir = 0; ir < neqn; ir++)
			{
				int indexi = neqn*N + ir;
				for(int i = 0; i < ndim; i++)
				{
					for (int beta = 0; beta < nnod; beta++)
						element->rhs( indexi) += w_i *p(beta, i, ir)*dPhi(N, i)*Phi(beta)*J;

					for(int j = 0; j < ndim; j++)
					{
						for (int beta = 0; beta < nnod; beta++)
							element->rhs( indexi) += w_i *q(beta, i, j, ir)*dPhi(N, i)*dPhi(beta,j)*J;
					}

					// += s(beta, r)
				}

				for(int M = 0; M < nnod; M++)
				{	
					for (int is = 0; is < neqn; is++)
					{
						int indexj = neqn*M+is;

						element->R( indexi, indexj) += w_i *eta(ir, is)*Phi(N)*Phi(M) * J;

						for (int i = 0; i < ndim; i++)
						{
							element->R( indexi, indexj)  += w_i * g(i, ir, is)* dPhi(N, i)* Phi(M)* J;

							for(int j=0; j<ndim; j++)
							{
								element->R( indexi, indexj) += w_i * h(i, j, ir, is) * dPhi(N, i) * dPhi(M, j)* J;
								element->R( indexi, indexj) += w_i * nd(i, ir, is)* Psi(N) * dPhi(M, j) * J;	// numerical diffusion
							}
						}						
					}
				}
			}
		}
	}


		/// Evaluate boundary integrals, if present
	void boundary(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
		// Loop over each face.
		for (unsigned int iface = 0; iface < element->face.size(); iface++)			// loop over each face
		{
			if (element->face[iface]->bc != 0 && element->face[iface]->bc < 100)
			{
				evalboundary(/*element->face[iface],*/ w, zeta, element->face[iface]->bc, element->face[iface]->n);
			}
		}
	}

	void evalboundary(/*Face * face,*/ Tensor<double, 1> & w, Tensor<double, 1> & zeta, int & ibndcnd, Tensor<int, 1> & n)
	{
		std::vector<Node *> node;
		for (int i=0; i<nbnod; i++)
		{
			node.push_back( element->node[n(i)] );
		}

		// Test Functions, normals
		Tensor<double, 1> phibou(nbnod);	// nbnod?
		Tensor<double, 1> en(ndim);
		Tensor<double, 2> dpidx(nbnod,ndim);

		if(ndim == 2)
		{
			// No changes necessary for adaptive meshing.
			TFNorm2D(node, phibou, en, dpidx, w, zeta);
		}
		else if (ndim == 3)
		{
			TFNorm3D(node, phibou, en, dpidx, w, zeta, element);

			// transform phibou, dpxi if adaptive
			// should be transforming dPhi/deta, not dPhi/dx. This is incorrect.
			// -> dPhi/dx is fixed already. Phibou is not. So don't fix dpidx.
			// in the future, need to fix phi, dphi in TFNorm3D or something to make it more clear.
			if (element->adap)	
			{
				Tensor<double, 2> fake(4,3);
				adap_face(phibou, fake, node, element);
			}
		}

		// Calculate values at gaussian point, along with k and mu
		sum(node,phibou,dpidx);

		// apply boundary conditions
		// Calculates tau and applies BC to tau, if any.
		if (ndim == 2)
		{
			eval1dboundary(ibndcnd, flow, node, en);
		}
		else if (ndim == 3)
		{
			eval2dboundary(ibndcnd, flow, en);
		}
		else if (ndim == 1)
		{
			flow.dT.clear();
			flow.dU.clear();
			flow.dV.clear();
			flow.tau.clear();
		}

		// Calculate a, bi, cij jacobians
		NS->calc_a(flow, a);
		NS->calc_b(flow, b);
		NS->calc_c(flow, c);
		NS->calc_d(d);

		sumFGB(node, phibou, dpidx);

		// Sum terms and calculate LHS, RHS
		sumBoundaryTerms();

		sumBoundaryMatrix(phibou, en, dpidx, n, ibndcnd);
	}

		/// Sum domain intermediate matrices
	void sumDomainTerms()
	{
		g.clear();					
		h.clear();
		p.clear();
		q.clear();
		s.clear();
		eta.clear();

		nd.clear();

		Tensor<double, 2> delta = Kroneker_delta<double,5>();

		for (int beta = 0; beta < nnod; beta++)
		{
			for(int i = 0; i < ndim; i++)
			{
				for( int r = 0; r < neqn; r++)
				{
					p(beta, i, r) += dt * ((*F[beta])(i, r) + (*G[beta])(i, r));

					for (int s=0; s<neqn; s++)
					{
						p(beta, i, r) += dtsq_2*(d(r, s) * ((*F[beta])(i, s) + (*G[beta])(i, s))) + dtsq_2*((a(i, r, s) + b(i, r, s)) * (*B[beta])(s));
					}
				}
			}
		}

		// g , nd
		for(int i = 0; i < ndim; i++)
		{
			for(int r = 0; r < neqn; r++)
			{
				for(int s = 0; s < neqn; s++)
				{
					g(i, r, s)  -= dt * (s1 * a(i, r, s) + s3 * b(i, r, s));
					nd(i, r, s) += dt * (s1 * a(i, r, s) + s3 * b(i, r, s));
					for (int t=0; t<neqn; t++)
					{
						g(i, r, s) -= dtsq_2*(s2*(d(r, t) * a(i, t, s)) + s6*(d(r, t) * (a(i, t, s) + b(i, t, s))) + s4*(d(r, t) * b(i, t, s)));
					}
				}
			}
		}
		for (int r = 0; r < neqn; r++)
		{
			for (int s=0; s < neqn; s++)
			{
				eta(r, s) -= dt*s5*d(r, s);
				for (int t=0; t < neqn; t++)
				{
					eta(r, s) -= dtsq_2*s6*(d(r, t) * d(t, s));
				}
				eta(r,s) += delta(r, s);
			}
		}

		for (int beta = 0; beta < nnod; beta++)
		{
			for (int i = 0; i < ndim; i++)
			{
				for(int j = 0; j < ndim; j++)
				{
					for( int r = 0; r < neqn; r++)
					{
						for(int s = 0; s < neqn; s++)
						{
							q(beta, i, j, r) -= dtsq_2*(a(i, r, s) + b(i, r, s))*((*F[beta])(j, s) + (*G[beta])(j, s));
						}
					}
				}
			}
		}

		for(int i = 0; i < ndim; i++)
		{
			for (int j = 0; j < ndim; j++)
			{
				for(int r = 0; r < neqn; r++)
				{
					for(int s = 0; s < neqn; s++)
					{
						h(i, j, r, s) -= (dt * s3 * c(i, j, r, s));
						for(int t=0; t<neqn; t++)
						{
							h(i, j, r, s) -= (- dtsq_2*(s2*(a(i, r, t)*a(j, t, s) + b(i, r, t) * a(j, t, s)) 
							+ s4*((a(i, r, t)*b(j, t, s) + b(i, r, t) * b(j, t, s))/* - (d(r, t) * c(i, j, t, s))*/) ));
						}
					}
				}
			}
		}

		for (int beta = 0; beta < nnod; beta++)
		{
			for (int r=0; r<neqn; r++)
			{
				s(beta, r) = dt * (*B[beta])(r);
				for (int ss=0; ss<neqn; ss++)
				{
					s(beta, r) += dtsq_2*(d(r, ss) * (*B[beta])(ss));
				}
			}
		}
	}

	void sumBoundaryTerms()
	{
		t.clear();
		u.clear();
		k.clear();
		m.clear();

		for (int beta = 0; beta < nbnod; beta++)
		{
			for (int i = 0; i < ndim; i++)
			{
				for(int r = 0; r < neqn; r++)
				{
					t(beta, i, r) = - dt * ((*F[beta])(i, r) + (*G[beta])(i, r));
					for (int s = 0; s < neqn; s++)
					{
						t(beta, i, r) -=  dtsq_2*(d(r, s) * ( (*F[beta])(i, s) + (*G[beta])(i, s))) - dtsq_2*((a(i, r, s) + b(i, r, s)) * (*B[beta])(s));
					}
				}
			}
		}

		for (int beta = 0; beta < nbnod; beta++)
		{
			for (int j = 0; j < ndim; j++)
			{
				for(int r = 0; r < neqn; r++)
				{
					for (int i = 0; i < ndim; i++)
					{
						for (int s = 0; s < neqn; s++)
						{
							u(beta, j, r, i) += dtsq_2 * (a(i, r, s) + b(i, r, s)) * ((*F[beta])(j, s) + (*G[beta])(j, s));
						}
					}
				}
			}
		}

		for(int i = 0; i < ndim; i++)
		{
			for(int r = 0; r < neqn; r++)
			{
				for(int s = 0; s < neqn; s++)
				{
					k(i, r, s) = dt * (s1 * a(i, r, s) + s3 * b(i, r, s)); 
					for (int t=0; t<neqn; t++)
					{
						k(i, r, s) +=dtsq_2*( s2*(d(r, t) * a(i, t, s)) 
								+ s6*(d(r, t) * (a(i, t, s) + b(i, t, s))) + s4*(d(r, t) * b(i, t, s)) );
					}
				}
			}
		}

		for(int i = 0; i < ndim; i++)
		{
			for (int j = 0; j < ndim; j++)
			{
				for(int r = 0; r < neqn; r++)
				{
					for(int s = 0; s < neqn; s++)
					{
						m(i, j, r, s) += (dt * s3 * c(i, j, r, s));
						for(int t = 0; t < neqn; t++)
						{
							m(i, j, r, s) -= dtsq_2*(s2*(a(i, r, t)*a(j, t, s) + b(i, r, t) * a(j, t, s)) 
								+ s4*((a(i, r, t)*b(j, t, s)) + (b(i, r, t) * b(j, t, s)) /*- (d(r, t) * c(i, j, t, s))*/));
						}
					}
				}
			}
		}
	}

	void sumBoundaryMatrix(Tensor<double, 1> & phibou, Tensor<double, 1> & en, Tensor<double, 2> & dpidx, Tensor<int, 1> & n, int ibndcnd)
	{

			// RHS
		for (int N = 0; N < nbnod; N++)
		{
			for (int ir = 0; ir < neqn; ir++)
			{
				int indexi = neqn*n(N)+ir;
				for(int i = 0; i < ndim; i++)
				{
					for (int beta = 0; beta < nbnod; beta++)
						element->rhs(indexi) += t(beta, i, ir) * en(i) * phibou(N)*phibou(beta);

					for(int j = 0; j < ndim; j++)
					{
						for (int beta = 0; beta < nbnod; beta++)
							element->rhs(indexi) += u(beta, j, ir, i) * en(i) * phibou(N)*dpidx(beta,j);
					}
				}

				for (int M = 0; M < nbnod; M++)
				{
					for( int is = 0; is < neqn; is++)
					{
						int indexj = neqn*n(M)+is;
						for (int i = 0; i < ndim; i++)
						{
							element->R(indexi, indexj) += k(i, ir, is) * phibou(N) * phibou(M) * en(i);
							for (int j = 0; j < ndim; j++)
							{
								element->R(indexi, indexj) += m(i, j, ir, is) * dpidx(M, j) * phibou(N) * en(i);
							}
						}
					}
				}
			}
		}
	}

		/// Calculates variables at the Gaussian point of interest, given the points
		/// on the line/area/volume and the appropriate test functions. 
		/// Therefore, it is universal for 1D, 2D and 3D.
	void sum(std::vector<Node *> & node, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi)
	{
		Thermo & thermo = Thermo::Instance();

		flow.clear();

		for (unsigned int N = 0; N < node.size(); N++)
		{
			flow.E   += Phi(N) * node[N]->E;
			flow.T   += Phi(N) * node[N]->T;

			// U, dU
			for (int i = 0; i < neqn; i++)
			{
				flow.U(i) += Phi(N) * node[N]->U(i);
				for (int j = 0; j < ndim; j++)
				{
					flow.dU(i, j) += dPhi(N, j) * node[N]->U(i);
				}
			}

			for (int i = 0; i < ndim; i++)
			{
				flow.x(i)  += Phi(N)     * node[N]->x(i);
				flow.v(i)  += Phi(N)     * node[N]->v(i);
				flow.dT(i) += dPhi(N, i) * node[N]->T;

				for (int dx1 = 0; dx1 < ndim; dx1++)
				{
					flow.dV(i, dx1) += dPhi(N, dx1) * node[N]->v(i);
				}
			}
		}

		flow.mu = thermo.calc_mu( flow.T);
		flow.k  = thermo.calc_k( flow.mu);
	}

		// we sum after BC's. Check about GWH. This may null some of his 
		// stuff where he rebuilds F, G, and zero's specific F, G terms.
	void sumFGB(std::vector<Node *> & node, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi)
	{
		Thermo & thermo = Thermo::Instance();

		for (unsigned int i=0; i < node.size(); i++)
		{
			NS->calc_F(node[i]->U, node[i]->p, *F[i]);
			NS->calc_G(flow.tau, node[i]->v, flow.dT, flow.k, thermo.creyn, *G[i]);
			NS->calc_B(node[i], *B[i]);
		}
	}

	T element;							/// Element pointer
	Flow flow;							/// Flow data structure
	NavierStokes * NS;

	GaussQuad<FDVGalerkin> gq;			/// 1D Gaussian Quadrature for boundary integral (2D)

	Tensor<double, 3> a;						/// dF/dU
	Tensor<double, 3> b;						/// dG/dU
	Tensor<double, 4> c;						/// dG/dU,j
	Tensor<double, 2> d;						/// dB/dU

	std::vector<Tensor<double, 1> *> B;		/// Source Term
	std::vector<Tensor<double, 2> *> F;		/// Convective Flux
	std::vector<Tensor<double, 2> *> G;		/// Viscous Flux

		/// domain intermediate matrices
	Tensor<double,2> s, eta;
	Tensor<double,3> p, g;
	Tensor<double,4> q, h;
	Tensor<double,3> nd;					/// numerical diffusion integral (GPG)

		/// boundary intermediate matrices
	Tensor<double,3> t, k;
	Tensor<double,4> u, m;

	double dt;
	double dtsq_2;

		/// dimensions
	const int nnod;
	const int nbnod;
	const int neqn;
	const int ndim;

		/// flowfield - dependant parameters
	double s1;
	double s2;
	double s3;
	double s4;
	double s5;
	double s6;
};

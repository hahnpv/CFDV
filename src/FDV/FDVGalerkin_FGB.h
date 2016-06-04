#pragma once

#include <iostream>
#include <numeric>
using namespace std;

	// 2D
#include "FDV/NavierStokes2D.h"
#include "FDV/Function/eval1dboundary.h"
#include "FiniteElement/TestFunctions1D.h"

	// 3D
#include "FDV/NavierStokes3D.h"
#include "FiniteElement/calc_TF_norm3D.h"					// norm + test functions
#include "FDV/Function/eval2dboundary.h"

	// Triangle FE
#include "FiniteElement/TestFunctionsTriangle2D.h"

	// Quadrature
#include "FiniteElement/GaussQuad.h"

#include "FDV/Function/FDVParam.h"
	// Adaptive Meshing
#include "adap/MeshRefine2D.h"

#include "FDV/Flow.h"

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

template<class T, class NavierStokes >
struct FDVGalerkin : public unary_function<T, void>
{
	FDVGalerkin(const int nnod, const int neqn, const int ndim, const int nbnod, double dt)
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
			nbnod(nbnod),				// number of boundary nodes, = 1/2 the number of nodes in 2D and 3D.
			gq(2),	
			dt(dt),
			dtsq_2(pow(dt,2.0)/2.0),
			a(ndim,neqn,neqn),
			b(ndim,neqn,neqn),
			c(ndim,ndim,neqn,neqn),
			d(neqn,neqn)
	{
		B.resize(nnod);
		G.resize(nnod);		// try Matrix1<G_NS>?
		F.resize(nnod);

		for (int i=0; i<nnod; i++)
		{
			F[i] = new Tensor<double, 2>(ndim, neqn);
			G[i] = new Tensor<double, 2>(ndim, neqn);
			B[i] = new Tensor<double, 1>(neqn);
		}
	}
	~FDVGalerkin()
	{
		delete NS;

		for (int i=0; i<nnod; i++)
		{
			delete F[i];
			delete G[i];
			delete B[i];
		}		
	}

	FDVGalerkin<T, NavierStokes> & operator[](double deltat)
	{
		dt = deltat;
		dtsq_2 = pow(dt, 2.0) / 2.0;
		return *this;
	}

	typedef void result_type;					// result_type for tr1 wrapper.
	result_type operator() (T e) 
	{
		element = e;

		FDVParam<Element *> fdvparam;						/// Calculate the FDV parameters

		fdvparam(e, s1, s2, s3, s4);
		s5 = 0.0;
		s6 = 0.0;

		// 2D Element Quadrature
		if (ndim == 2)
		{
			gq.two(this,&FDVGalerkin<T, NavierStokes>::domain);
			gq.one(this,&FDVGalerkin<T, NavierStokes>::boundary);
		}

		// 3D Element Quadrature
		else if (ndim == 3)
		{
			gq.three(this,&FDVGalerkin<T, NavierStokes>::domain);
			gq.two(this,&FDVGalerkin<T, NavierStokes>::boundary);
		}
	}

		/// Calculate the domain integral
	void domain(Tensor<double, 1> & w, Tensor<double, 1> & zeta)
	{
		std::vector<Node *> & node = element->node;

//		PrecisionTimer timr;
//		timr.start();

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
		else if (ndim == 2)
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

//		timr.stop();
//		fout << timr.read() << ", ";
//		timr.start();
		
		// 2. Calculate flow variables, viscocity, k, tau
		sum(node,Phi,dPhi);
		flow.tau.calculate(flow.mu,flow.dV);

//		timr.stop();
//		fout << timr.read() << ", ";
//		timr.start();
		
		// 4. Calculate a, b, c jacobians and F,G fluxes
		NS->calc_a(flow, a);		//		a.update(flow);
		NS->calc_b(flow, b);		//		b.update(flow);
		NS->calc_c(flow, c);		//		c.update(flow);
		NS->calc_d(d);

		sumFGB(node, Phi, dPhi);

//		timr.stop();
//		fout << timr.read() << ", ";
//		timr.start();
		
		// 5. Assemble vectors 
		sumDomainTerms();

//		timr.stop();
//		fout << timr.read() << ", ";
//		timr.start();
		
		sumDomainMatrix(w, Phi, dPhi, J);

//		timr.stop();
//		fout << timr.read() << endl;
	}

	void sumDomainMatrix(Tensor<double, 1> & w, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi, double & J)
	{
		double w_i = 1;
		for (unsigned int i=0; i<w.size(0); i++)
			w_i *= w(i);

		double twopiy = 1.0;
		if (true)
		{
		//	twopiy = flow.x(1);
		}

			// RHS
		for (int N = 0; N < nnod; N++)
		{
			for (int ir = 0; ir < neqn; ir++)
			{
				int rhs_N = neqn*N + ir;
				for(int i = 0; i < ndim; i++)
				{
					for(int beta = 0; beta < nnod; beta++)
					{
						element->rhs(rhs_N) += w_i *p(beta, i, ir)*Phi(beta)*dPhi(N, i)*J * twopiy;					// ok
				
						for(int j = 0; j < ndim; j++)
						{
							element->rhs(rhs_N) += w_i *q(beta,i, j, ir)*dPhi(N, i)*dPhi(beta, j)*J * twopiy;		// ok
						}
		
						element->rhs(rhs_N) += w_i * s(beta, ir) * Phi(N) * Phi(beta) * J * twopiy;
						// += s(beta, r)
					}
				}

				for(int M = 0; M < nnod; M++)
				{	
					for (int is = 0; is < neqn; is++)
					{
						int rhs_M = neqn*M+is;

						element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M) += w_i *eta(ir, is)*Phi(N)*Phi(M) * J * twopiy;

						for (int i = 0; i < ndim; i++)
						{
							element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M)  += w_i *g(i, ir, is)* dPhi(N, i)* Phi(M)* J * twopiy;

							for(int j=0; j<ndim; j++)
							{
								element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M) += w_i *h(i, j, ir, is) * dPhi(N, i) * dPhi(M, j)* J * twopiy;
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
				evalboundary(element->face[iface], w, zeta, element->face[iface]->bc, element->face[iface]->n);
			}
		}
	}

	void evalboundary(Face * face, Tensor<double, 1> & w, Tensor<double, 1> & zeta, int & ibndcnd, Tensor<int, 1> & n)
	{
		std::vector<Node *> node;
		for (int i=0; i<nbnod; i++)
		{
			node.push_back( element->node[n(i)] );
		}

//		PrecisionTimer timr;
//		timr.start();

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

//		timr.stop();
//		foutbc << timr.read() << ", ";
//		timr.start();

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

		// Calculate a, bi, cij jacobians
		NS->calc_a(flow, a);		//		a.update(flow);
		NS->calc_b(flow, b);		//		b.update(flow);
		NS->calc_c(flow, c);		//		c.update(flow);
		NS->calc_d(d);

//		sumFGB(node, phibou, dpidx);

		// Sum terms and calculate LHS, RHS
		sumBoundaryTerms();						// here

		sumBoundaryMatrix(phibou, en, dpidx, n, ibndcnd);

//		timr.stop();
//		foutbc << timr.read() << endl;

	}

		/// Sum domain intermediate matrices - OK except B term (new sum).
	void sumDomainTerms()
	{
		g.clear();					
		h.clear();
		p.clear();
		q.clear();
		s.clear();
		eta.clear();

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
			for(int r = 0; r < neqn; r++)
			{
				for(int s = 0; s < neqn; s++)
				{
					g(i, r, s) -= dt * (s1 * a(i, r, s) + s3 * b(i, r, s));
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
				if (r == s)				// Kroneker delta
				{
					eta(r, s) += 1.0;
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

			// not implemented yet; needs fix.
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
		double twopiy = 1.0;
		if (true)
		{
		//	twopiy =  flow.x(1);
		}
	 
			// RHS
		for (int N = 0; N < nbnod; N++)
		{
			for (int ir = 0; ir < neqn; ir++)
			{
				int indexi = neqn*n(N)+ir;
				for(int i = 0; i < ndim; i++)
				{
					for (int beta = 0; beta < nbnod; beta++)
					{
						element->rhs(indexi) += t(beta, i, ir) * en(i) * phibou(N) * phibou(beta) * twopiy;
						
						for(int j = 0; j < ndim; j++)
						{
							element->rhs(indexi) += u(beta, j, ir, i) * en(i) * phibou(N) * dpidx(beta,j) * twopiy;
						}
					}
				}

				for (int M = 0; M < nbnod; M++)
				{
					for( int is = 0; is < neqn; is++)
					{
						int indexj = neqn*n(M)+is;
						for (int i = 0; i < ndim; i++)
						{
							element->R(indexi, indexj) += k(i, ir, is) * phibou(N) * phibou(M) * en(i) * twopiy;
							for (int j = 0; j < ndim; j++)
							{
								element->R(indexi, indexj) += m(i, j, ir, is) * dpidx(M, j) * phibou(N) * en(i) * twopiy;
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
		// B should go here eventually.
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
	Tensor<double,2> eta;
	Tensor<double,3> g, p;
	Tensor<double,4> q, h;

		/// boundary intermediate matrices
	Tensor<double,2> s;
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

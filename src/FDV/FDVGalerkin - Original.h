#pragma once

#include <iostream>
#include <numeric>
using namespace std;

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
	// 3. Someday: don't need to pre-sum F, G, can do in-line...

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
			p(ndim,neqn),				// beta already summed
			q(ndim,ndim,neqn),			// beta already summed
			eta(neqn,neqn),
			t(ndim,neqn),				// beta already summed
			u(ndim,neqn,ndim),
			k(ndim,neqn,neqn),
			m(ndim,ndim,neqn,neqn),
			s(neqn),					// beta already summed
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
		G.resize(nnod);
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

		FDVParam<Element *> fdvparam;			/// Calculate the FDV parameters

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
		
		// 2. Calculate flow variables, viscocity, k, tau
		sum(node,Phi,dPhi);
		flow.tau.calculate(flow.mu,flow.dV);

		// 4. Calculate a, b, c jacobians and F,G fluxes
		NS->calc_a(flow, a);		//		a.update(flow);
		NS->calc_b(flow, b);		//		b.update(flow);
		NS->calc_c(flow, c);		//		c.update(flow);
		NS->calc_d(d);

		sumFG(node, Phi, dPhi);

		// 5. Assemble vectors 
		sumDomainTerms();

		sumDomainMatrix(w, Phi, dPhi, J);
	}

	void sumDomainMatrix(Tensor<double, 1> & w, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi, double & J)
	{
		double w_i = 1;
		for (unsigned int i=0; i<w.imax(); i++)
			w_i *= w(i);

		// i think the problem is here, if I revert to "correct" method here, but incorrect method
		// in bc, breaks.
		double twopiy = 1.0;
//		if (true)
//		{
		//	twopiy = /*2.0*3.14159**/flow.x(1);
//		}


			// RHS
		for (int N = 0; N < nnod; N++)
		{
		//	twopiy = element->node[N]->x(1);
			for (int ir = 0; ir < neqn; ir++)
			{
				int rhs_N = neqn*N + ir;
				for(int i = 0; i < ndim; i++)
				{
					element->rhs(/*neqn*N+ir*/ rhs_N) += w_i *p(i, ir)*dPhi(N, i)*J*twopiy;

					for(int j = 0; j < ndim; j++)
					{
						element->rhs(/*neqn*N+ir*/ rhs_N) += w_i *q(i, j, ir)*dPhi(N, i)*J*twopiy;
					}

					// += s(beta, r)
				}
//			}				// Making N, ir one loop results in a 2% speedup overall in FDV (7% in function)
//		}

		// LHS
//		for(int N = 0; N < nnod; N++)
//		{
//			for (int ir = 0; ir < neqn; ir++)
//			{
//				int rhs_N = neqn*N + ir;
				for(int M = 0; M < nnod; M++)
				{	
					for (int is = 0; is < neqn; is++)
					{
						int rhs_M = neqn*M+is;

						element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M) += w_i *eta(ir, is)*Phi(N)*Phi(M) * J*twopiy;

						for (int i = 0; i < ndim; i++)
						{
							element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M)  += w_i *g(i, ir, is)* dPhi(N, i)* Phi(M)* J*twopiy;

							for(int j=0; j<ndim; j++)
							{
								element->R(/*neqn*N+ir*/ rhs_N, /*neqn*M+is*/ rhs_M) += w_i *h(i, j, ir, is) * dPhi(N, i) * dPhi(M, j)* J*twopiy;
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

		// Calculate a, bi, cij jacobians
		NS->calc_a(flow, a);		//		a.update(flow);
		NS->calc_b(flow, b);		//		b.update(flow);
		NS->calc_c(flow, c);		//		c.update(flow);
		NS->calc_d(d);

		sumFG(node, phibou, dpidx);

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

		Tensor<double, 2> delta = Kroneker_delta<double,5>();

		for(int i = 0; i < ndim; i++)
		{
			for( int r = 0; r < neqn; r++)
			{
				p(i, r) += dt * (flow.F(i, r) + flow.G(i, r));
				// outside of loop
				for (int s=0; s<neqn; s++)
				{
					p(i, r) += dtsq_2*(d(r, s) * (flow.F(i, s) + flow.G(i, s))) + dtsq_2*((a(i, r, s) + b(i, r, s)) * flow.B(s));
				}
			}
		}

		// g 
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
/*				if (r == s)				// Kroneker delta
				{
					eta(r, s) += 1.0;
				}
*/
				eta(r,s) += delta(r, s);
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
						q(i, j, r) -= dtsq_2*(a(i, r, s) + b(i, r, s))*(flow.dF(j, s) + flow.dG(j, s));
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

		for (int r=0; r<neqn; r++)
		{
			s(r) = dt * flow.B(r);
			for (int ss=0; ss<neqn; ss++)
			{
				s(r) += dtsq_2*(d(r, ss) * flow.B(ss));
			}
		}
	}

	void sumBoundaryTerms()
	{
		t.clear();
		u.clear();
		k.clear();
		m.clear();

		for (int i = 0; i < ndim; i++)
		{
			for(int r = 0; r < neqn; r++)
			{
				t(i, r) = - dt * (flow.F(i, r) + flow.G(i, r));
				for (int s = 0; s < neqn; s++)
				{
					t(i, r) -=  dtsq_2*(d(r, s) * ( flow.F(i, s) + flow.G(i, s))) - dtsq_2*((a(i, r, s) + b(i, r, s)) * flow.B(s));
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
						u(j, r, i) += dtsq_2 * (a(i, r, s) + b(i, r, s)) * (flow.dF(j, s) + flow.dG(j, s));
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
//		if (true)
//		{
		//	twopiy = /*2.0*3.14159**/flow.x(1);
//		}

		// can try summing x here, should be identical to Flow ... 

			// RHS
		for (int N = 0; N < nbnod; N++)
		{
		//	twopiy = element->node[n(N)]->x(1);
			for (int ir = 0; ir < neqn; ir++)
			{
				int indexi = neqn*n(N)+ir;
				for(int i = 0; i < ndim; i++)
				{
					element->rhs(indexi) += t(i, ir) * en(i) * phibou(N)*twopiy;

					for(int j = 0; j < ndim; j++)
					{
						element->rhs(indexi) += u(j, ir, i) * en(i) * phibou(N)*twopiy;
					}
				}
//			}
//		}


//		for (int N = 0; N < nbnod; N++)
//		{
//			for (int ir = 0; ir < neqn; ir++)
//			{
				for (int M = 0; M < nbnod; M++)
				{
					for( int is = 0; is < neqn; is++)
					{
						int indexj = neqn*n(M)+is;
						for (int i = 0; i < ndim; i++)
						{
							element->R(indexi, indexj) += k(i, ir, is) * phibou(N) * phibou(M) * en(i)*twopiy;
							for (int j = 0; j < ndim; j++)
							{
								element->R(indexi, indexj) += m(i, j, ir, is) * dpidx(M, j) * phibou(N) * en(i)*twopiy;
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
	void sumFG(std::vector<Node *> & node, Tensor<double, 1> & Phi, Tensor<double, 2> & dPhi)
	{
		Thermo & thermo = Thermo::Instance();

		for (unsigned int i=0; i < node.size(); i++)
		{
			NS->calc_F(node[i]->U, node[i]->p, *F[i]);
			NS->calc_G(flow.tau, node[i]->v, flow.dT, flow.k, thermo.creyn, *G[i]);
			NS->calc_B(node[i], *B[i]);
		}
			// then calc flow.f, dF here
		for (unsigned int N = 0; N < node.size(); N++)
		{
			for (int i = 0; i < neqn; i++)
			{
				for (int j = 0; j < ndim; j++)
				{
					flow.G(j, i)  += Phi(N)     * (*G[N])(j, i);
					flow.dG(j, i) += dPhi(N, j) * (*G[N])(j, i);
					flow.F(j, i)  += Phi(N) *     (*F[N])(j, i);
					flow.dF(j, i) += dPhi(N, j) * (*F[N])(j, i);
					flow.B(j)/*[i]*/  += Phi(N) * (*B[N])(j)/*[i]*/;		// b subscripts are bad ... 
					flow.dB(j, i) += dPhi(N, j) * (*B[N])(j)/*[i]*/;
				}
			}
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
	Tensor<double,2> p, eta;
	Tensor<double,3> g, q;
	Tensor<double,4> h;

		/// boundary intermediate matrices
	Tensor<double,1> s;
	Tensor<double,2> t;
	Tensor<double,3> u, k;
	Tensor<double,4> m;

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

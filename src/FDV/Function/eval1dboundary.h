#pragma once

//	Function to evaluate 1D boundary for 2D CFD
//	Need to re-visit all of the cases here, I think they are all suspect.

#include "../Flow.h"

// phibou not used ... 
void eval1dboundary(int ibndcnd, Flow & flow, std::vector<Node *> & node, Tensor<double, 1> & en)		// nodes temp test
{
		Tensor<double, 1> Vdotn(2);				// V dot n
		Tensor<double, 1> n = unit(en);			// normal, where ||n|| = 1.
//		cout << "en / n: " << en << " " << n << endl;

		switch ( ibndcnd )						// Case 2 validated
		{
			case 2 :							// symmetry
			case 12:
/*				flow.v(1)		= 0.0;
				flow.dV(0, 1) = 0.0;
				flow.dV(1, 1) = 0.0;
				flow.dT(1)    = 0.0;
				for (unsigned int i=0; i<flow.dU.imax(); i++)
				{
					flow.dU(i, 1) = 0.0;		// my dU is reverse his, make sure to observe below in equation formulations.
				}
				//--> tau(0,1) == tau(1,0) == 0, applied after tau is calculated.
*/
				flow.v(0) -= dot(flow.v, n) * n(0);			// v_n = 0
				flow.v(1) -= dot(flow.v, n) * n(1);

				flow.dT(0) -= dot(flow.dT, n) * n(0);		// dT_n = 0
				flow.dT(1) -= dot(flow.dT, n) * n(1);

				// flow.dV(tangent)_dn = 0				-> need to figure out.
			
			break;

			case 1 :							// noslip wall, u,v=0
			case 51:							// T=Twall or dTdn=0 (constant wall temp or adiabatic)
			/*
				flow.v.clear();
				flow.dV(0, 1) = 0.0;
				flow.dV(1, 1) = 0.0;
				flow.dT(1)    = 0.0;
				for (unsigned int i=0; i<flow.dU.imax(); i++)			// number of equations
				{
					flow.dU(i, 1) = 0.0;
				}
			*/
				// FIXME adiabatic only?
				flow.dT(0) -= dot(flow.dT, n) * n(0);
				flow.dT(1) -= dot(flow.dT, n) * n(1);

				// flow.dV = 0 since v at all nodes = 0.
				// pressure gradient?	

			break;

			case -1:							// inlet/exit (also 6)
				{
					flow.dV(0, 0) = 0.0;
					flow.dV(1, 0) = 0.0;

					flow.dV(0, 1) = 0.0;
					flow.dV(1, 1) = 0.0;

					flow.dT(0)    = 0.0;
//					for (unsigned int i=0; i<flow.dU.imax(); i++)
//					{
//						flow.dU(i, 0)    = 0.0;		// my dU is reverse his, make sure to observe below in equation formulations.
//					}
				}
				break;
			case  6:							// 
				flow.dV(0, 0) = 0.0;
				flow.dV(1, 0) = 0.0;

				flow.dV(0, 1) = 0.0;			// TEST, subsonic
				flow.dV(1, 1) = 0.0;

				flow.dT(0)    = 0.0;

				for (unsigned int i=0; i<flow.dU.imax(); i++)
				{
					flow.dU(i, 0)    = 0.0;		// my dU is reverse his, make sure to observe below in equation formulations.
				}
			break;
			
			case -2:
			case  7:							// top surface [normal derivatives =0]
												// this is the same as reflection except no v[1]=0, and no tau changes.
				flow.dV(0, 1) = 0.0;
				flow.dV(1, 1) = 0.0;
				flow.dT(1)    = 0.0;
				for (unsigned int i=0; i<flow.dU.imax(); i++)			// # of equations
				{
					flow.dU(i, 1) = 0.0;		// my dU is reverse his, make sure to observe below in equation formulations.
				}
			break;

			case 4:								// no-slip compression corner (updated)
			{									// verify, converted to Tensor<> without a test case.
				Tensor<double, 1> n(2);			// try to condense dot's and such
				Tensor<double, 1> dUdot(4);
				Tensor<double, 1> dFdot(4);

				double dudot;
				double dvdot;
				double dtdot;

				flow.v.clear();

				n(0) = en(0);
				n(1) = en(1);

				// (AdotB)*B = projection onto B
				// B - projection = B with projection vector zero'd
				dudot = (flow.dV(0, 0) * n(0) + flow.dV(0, 1) * n(1));
				dvdot = (flow.dV(1, 0) * n(0) + flow.dV(1, 1) * n(1));
				dtdot = (flow.dT(0) * n(0) + flow.dT(1) * n(1));

				flow.dV(0, 0) -= dudot * n(0);
				flow.dV(0, 1) -= dudot * n(1);

				flow.dV(1, 0) -= dvdot * n(0);
				flow.dV(1, 1) -= dvdot * n(1);

				flow.dT(0) -= dtdot * n(0);
				flow.dT(1) -= dtdot * n(1);

				for (unsigned int i=0; i<flow.dU.imax(); i++)
				{
					dUdot(i) = flow.dU(i, 0) * n(0) + flow.dU(i, 1) * n(1);
		//			dFdot[i] = flow.dF[0][i] * n[0] + flow.dF[1][i] * n[1];
				}

					for (unsigned int i=0; i<flow.dU.imax(); i++)
					{
						for (int j = 0; j < 2; j++)
						{
							flow.dU(i, j) -= dUdot(i) * n(j);
			//				flow.dF[j][i] -= dFdot[i] * n[j];
						}
					}
			}
			break;

			case -12:			// top of lid driven cavity
				flow.dV(0, 1) = 0.0; 
				flow.dV(1, 1) = 0.0; 
				flow.dT(1)    = 0.0; 
				flow.v(1)    = 0.0; 
				for (unsigned int i=0; i<flow.dU.imax(); i++)		// # eqns
				{	
					flow.dU(i, 1) = 0.0;
				}		
			break;


			case 11:							// side walls of lid driven cavity
			case 21:
				flow.v.clear();
				flow.dV(0, 0) = 0.0;
				flow.dV(1, 0) = 0.0;
				flow.dT(1)    = 0.0;
				for(unsigned int i=0; i<flow.dU.imax(); i++)			// eqns
				{
					flow.dU(i, 1) = 0.0;
				}
			break;

			case 8:								// compression corner exit (derivs normal to flow =0)
			{									// unverified / condense dots / etc					
				Tensor<double, 1> n(2);
				Tensor<double, 1> dUdot(4);
				Tensor<double, 1> dFdot(4);

				double dudot;
				double dvdot;
				double dtdot;
				double dl;

				double alpha;
				alpha = 10 * (3.14159/180.0);		// 10 degree ramp angle
				
				dl= abs(node[1]->x(1) - node[0]->x(1));

				n(0) = cos(alpha) * dl;
				n(1) = sin(alpha) * dl;


				// (AdotB)*B = projection onto B
				// B - projection = B with projection vector zero'd
	
				dudot = (flow.dV(0, 0) * n(0) + flow.dV(0, 1) * n(1));
				dvdot = (flow.dV(1, 0) * n(0) + flow.dV(1, 1) * n(1));
				dtdot = (flow.dT(0) * n(0) + flow.dT(1) * n(1));

				flow.dV(0, 0) -= dudot * n(0);
				flow.dV(0, 1) -= dudot * n(1);

				flow.dV(1, 0) -= dvdot * n(0);
				flow.dV(1, 1) -= dvdot * n(1);

				flow.dT(0) -= dtdot * n(0);
				flow.dT(1) -= dtdot * n(1);

				for (unsigned int i=0; i<flow.dU.imax(); i++)
				{
					dUdot(i) = flow.dU(i, 0) * n(0) + flow.dU(i, 1) * n(1);
		//			dFdot[i] = flow.dF[0][i] * n[0] + flow.dF[1][i] * n[1];
				}

					for (unsigned int i=0; i<flow.dU.imax(); i++)
					{
						for (int j = 0; j < 2; j++)
						{
							flow.dU(i, j) -= dUdot(i) * n(j);
		//					flow.dF[j][i] -= dFdot[i] * n[j];
						}
//					    dfxdy(iequa)=0.0				?? I have not seen this term before. Or is it antiquated, zero dF[1]?
					}
			}
			break;

			case 0:								// no BC
				cout << "you shouldn't be here!" << endl;
			break;

			default:							// catch boundary conditions yet to be implemented
				cout << "No handler for boundary condition " << ibndcnd << endl;
				cin.get();
			break;
		}

		flow.tau.calculate(flow.mu,flow.dV);

		// If one of the boundary conditions prescribes tau_ij == tau_ji == 0,
		// then fix it here.
		if (ibndcnd == 2 || ibndcnd == 12 || ibndcnd == 7)		// 6 == exit, subsonic test 
		{
			flow.tau(0,1) = 0.0;
			flow.tau(1,0) = 0.0;
		}
		if (ibndcnd == 6)
		{
			flow.tau(0, 0) = 0.;
			flow.tau(1, 0) = 0.;
		}
};


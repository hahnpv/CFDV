#pragma once

#include "../Flow.h"
#include <cmath>
/*
#include "FDV/BC/BoundaryConditionInflow.h"
#include "FDV/BC/BoundaryConditionOutflow.h"
#include "FDV/BC/BoundaryConditionPlate.h"
#include "FDV/BC/BoundaryConditionDissipation.h"
#include "FDV/BC/BoundaryConditionReflection.h"
*/
	// Evaluation of a 2D surface integration of a boundary face in 3D
	// Enough are implemented to do flat plate studies, need to implement lid flow and ramp ones as well.
void eval2dboundary(int ibndcnd, Flow & flow, Tensor<double, 1> & en)
{
	switch ( ibndcnd )
	{
		case 2 :							// symmetry
		case 12:							// Symmetry in x-y plane, w=tau13=tau23=q3=0 
		{
//			int neqn = flow.dU.imax;
//			BoundaryConditionReflection reflection(neqn, 2);
//			reflection.ApplyNeumann(flow, en);

//                        flow.v(2) = 0; dirichlet enforced.

                        for (int i = 0; i < 3; i++)
                        {
                                flow.dV(i, 2) = 0;
                        }
                        flow.dT(2) = 0;
		}
		break;

		case 22:						// symmetry in the X/Z plane NOTE DIFFERENT THAN gwh
		{								// Symmetry in y-z plane, u=tau12=tau13=q1=0 
//			int neqn = flow.dU.imax;
//			BoundaryConditionReflection reflection(neqn, 1);
//			reflection.ApplyNeumann(flow, en);
			for (int i = 0; i < 3; i++)
			{
				flow.dV(i, 1) = 0;
//				flow.dV(0, 1) = 0;
//				flow.dV(2, 1) = 0;
			}
			flow.dT(1) = 0;
			flow.tau.calculate( flow.mu, flow.dV);

			flow.tau(0,1) = 0.0;
			flow.tau(1,1) = 0.0;		// new
			flow.tau(2,1) = 0.0;
			flow.tau(1,0) = 0.0;
			flow.tau(1,1) = 0.0;		// new
			flow.tau(1,2) = 0.0;
		}
		break;

		case 1:
		case 51:
		{												//  No slip wall, u=v=w=0 
			if (Thermo::Instance().adiabatic)
				flow.dT(2) = 0.0;

			flow.tau.calculate(flow.mu, flow.dV);
		}
		break;

		case 21:				// nonslip wall XZ plane
		case 31:
		{
			if (Thermo::Instance().adiabatic)
				flow.dT(1) = 0.0;

			flow.tau.calculate(flow.mu, flow.dV);
		}
		break;

		case 11:		// Noslip Wall in the YZ plane
		{
			if (Thermo::Instance().adiabatic)
				flow.dT(0) = 0.0;

			flow.tau.calculate(flow.mu, flow.dV);
		}
		break;

		case -1:
		{
		//	int neqn = flow.dU.imax();
		//	BoundaryConditionInflow inflow(neqn);
		//	inflow.ApplyNeumann(flow, en);
		}
		return;

		case 6:
		{												// Inlet and Exit, normal derivatives = 0 
		//	int neqn = flow.dU.imax;
		//	BoundaryConditionOutflow outflow(neqn);
		//	outflow.ApplyNeumann(flow, en);

//                        flow.dV(0, 0) = 0.0;
  //                      flow.dV(1, 0) = 0.0;
    //                    flow.dV(2, 0) = 0.0;
                        flow.dT(0) = 0.0;
                        flow.tau.calculate(flow.mu, flow.dV);
			for (int i = 0; i < 3; i++)
			{
				flow.tau(i, 0) = 0;			// not valid if curved.
			}
		}
		break;

		case 7:								//  Normal derivatives on top surface = 0 
		{
		//	int neqn = flow.dU.imax;
		//	BoundaryConditionDissipation dissipation(neqn);
		//	dissipation.ApplyNeumann(flow, en);

                        flow.dT(2) = 0.0;
                        flow.tau.calculate(flow.mu, flow.dV);
                        for (int i = 0; i < 3; i++)
                        {
                                flow.tau(i, 2) = 0;                     // not valid if curved.
                        }

		}

		case 77:								//  Normal derivatives on side surface = 0 
		case 32:								// like 7 in the YZ plane
		{
//			for (int i = 0; i < 3; i++)
//			{
//				flow.dV(i, 1) = 0;
//			}
//			flow.dT(1) = 0;

                        flow.dT(1) = 0.0;
                        flow.tau.calculate(flow.mu, flow.dV);
                        for (int i = 0; i < 3; i++)
                        {
                                flow.tau(i, 1) = 0;                     // not valid if curved.
                        }

		}
		break;

		// should be converted to 3D, verify and use.
		// convert the hand products into loop, verify sameness.
		// optimally work out vector functions or whatever to generate normals. IE, dV dot normal instead of loop.
		case 4:								// no-slip compression corner (updated)
		{									// not verified, no available test case, condense dots, etc.
			Tensor<double, 1> n(3);			// re-name normal vec in header no need to copy it ...
			Tensor<double, 1> dUdot(5);

			double dudot;
			double dvdot;
			double dtdot;

			flow.v.clear();

			n(0) = en(0);
			n(1) = en(1);
			n(2) = en(2);

			// (AdotB)*B = projection onto B
			// B - projection = B with projection vector zero'd
			dudot = (flow.dV(0, 0) * n(0) + flow.dV(0, 1) * n(1) + flow.dV(0, 2) * n(2));
			dvdot = (flow.dV(1, 0) * n(0) + flow.dV(1, 1) * n(1) + flow.dV(1, 2) * n(2));
			dtdot = (flow.dT(0) * n(0) + flow.dT(1) * n(1) + flow.dT(2) * n(2));

			flow.dV(0, 0) -= dudot * n(0);
			flow.dV(0, 1) -= dudot * n(1);
			flow.dV(0, 2) -= dudot * n(2);

			flow.dV(1, 0) -= dvdot * n(0);
			flow.dV(1, 1) -= dvdot * n(1);
			flow.dV(1, 2) -= dvdot * n(2);

			flow.dT(0) -= dtdot * n(0);
			flow.dT(1) -= dtdot * n(1);
			flow.dT(2) -= dtdot * n(2);

			for (unsigned int i=0; i<flow.dU.imax(); i++)
			{
				dUdot(i) = flow.dU(i, 0) * n(0) + flow.dU(i, 1) * n(1) + flow.dU(i, 2) * n(2);
			}

			for (unsigned int i = 0; i < flow.dU.jmax(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					flow.dU(i, j) -= dUdot(i) * n(j);
				}
			}
		}
		break;

		// Need to convert to 3D
		// Omitted for now ...
		case 8:								// compression corner exit (derivs normal to flow =0)
		{									// SPECIFIC TO compression corner in XY plane (ie XY plate up in +/-Z)
//			flow.tau.calculate(flow.mu, flow.dV);
/*

			Tensor<double, 1> nv = flow.v;

			double nvmag, nmag;
			for(int i = 0; i < 3; i++)
			{
				nmag += en(i) * en(i);
				nvmag += nv(i) * nv(i);
			}
			nmag = sqrt(nmag);
			nvmag = sqrt(nvmag);

			double angle = acos( dot(nv, en) / (nmag*nvmag) ) * (180.0/3.14159);
*/
//			cout << "angle: " << angle << " en: " << en << " nv: " << nv << " mag: " << nmag << " " << nvmag << endl;
			// angle appears valid

/*
			double dTm = sqrt( flow.dT(0) * flow.dT(0) + flow.dT(1) * flow.dT(1) );

			flow.dT(0) -= dTm * sin(angle);
			flow.dT(1) -= dTm * cos(angle);


			double dVx = sqrt( flow.dV(0, 0) * flow.dV(0) + flow.dV(0, 1) * flow.dV(1) );

			flow.dV(0, 0) -= dVx * sin(angle);
			flow.dV(0, 1) -= dVx * cos(angle);

			double dVy = sqrt( flow.dV(1, 0) * flow.dV(0) + flow.dV(1, 1) * flow.dV(1) );

			flow.dV(1, 0) -= dVy * sin(angle);
			flow.dV(1, 1) -= dVy * cos(angle);
*/

/*
			for (unsigned int i=0; i<flow.dU.imax(); i++)
			{
				dUdot(i) = flow.dU(i, 0) * n(0) + flow.dU(i, 1) * n(1) + flow.dU(i, 2) * n(2);
			}

			for (unsigned int i = 0; i < flow.dU.size(i); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					flow.dU(i, j) -= dUdot(i) * n(j);
				}
			}
*/
/*
			int neqn = flow.dU.imax;
			BoundaryConditionOutflow outflow(neqn);
			outflow.ApplyNeumann(flow, en);
		*/
/*
			Tensor<double, 1> c(3);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					c(i) = flow.dV(i, j) * en(j);

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					flow.dV(i, j) -= c(i);

			flow.tau.calculate(flow.mu, flow.dV);
*/
/*
			for (int i = 0; i < 3; i++)
			{
				flow.dV(i, 0) = 0;
			}
			flow.dT(0) = 0;
			flow.tau.calculate(flow.mu, flow.dV);
*/
/*
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if (en(j) != 0)
						flow.tau(i, j) = 0;
//					flow.tau(i, j) -= flow.tau(i ,j) * en(j);	
				}
					if (en(i) != 0)
						flow.dT(i) = 0;

//				flow.dT(i) -= flow.dT(i) * en(i);
			}
			flow.tau.calculate(flow.mu, flow.dV);
*/
/*
				Tensor<double, 1> n(3);
				Tensor<double, 1> dUdot(5);
				Tensor<double, 1> dFdot(5);

				double dudot;
				double dvdot;
				double dwdot;
				double dtdot;

				n(0) = en(0);
				n(1) = en(1);
				n(2) = en(2);

				// (AdotB)*B = projection onto B
				// B - projection = B with projection vector zero'd
	
				dudot = (flow.dV(0, 0) * n(0) + flow.dV(0, 1) * n(1) + flow.dV(0, 2) * n(2));
				dvdot = (flow.dV(1, 0) * n(0) + flow.dV(1, 1) * n(1) + flow.dV(1, 2) * n(2));
				dwdot = (flow.dV(2, 0) * n(0) + flow.dV(2, 1) * n(1) + flow.dV(2, 2) * n(2));
				dtdot = (flow.dT(0) * n(0) + flow.dT(1) * n(1) + flow.dT(2) * n(2));

				flow.dV(0, 0) = dudot * n(0);
				flow.dV(0, 1) = dudot * n(1);
				flow.dV(0, 2) = dudot * n(2);

				flow.dV(1, 0) = dvdot * n(0);
				flow.dV(1, 1) = dvdot * n(1);
				flow.dV(1, 2) = dvdot * n(2);

				flow.dV(2, 0) = dwdot * n(0);
				flow.dV(2, 1) = dwdot * n(1);
				flow.dV(2, 2) = dwdot * n(2);

				flow.dT(0) = dtdot * n(0);
				flow.dT(1) = dtdot * n(1);
				flow.dT(2) = dtdot * n(1);
				*/
/*
				for (unsigned int i=0; i<flow.dU.imax(); i++)
				{
					dUdot(i) = flow.dU(i, 0) * n(0) + flow.dU(i, 1) * n(1) + flow.dU(i, 2) * n(2);
		//			dFdot[i] = flow.dF[0][i] * n[0] + flow.dF[1][i] * n[1];
				}

					for (unsigned int i=0; i<flow.dU.imax(); i++)
					{
						for (int j = 0; j < 3; j++)
						{
							flow.dU(i, j) -= dUdot(i) * n(j);
		//					flow.dF[j][i] -= dFdot[i] * n[j];
						}
//					    dfxdy(iequa)=0.0				?? I have not seen this term before. Or is it antiquated, zero dF[1]?
					}
//			}
*/

				flow.tau.calculate(flow.mu, flow.dV);

		}			// does the same to dF but if we construct with these product, should be OK?
		break;

		case -12:			// top of lid driven cavity, may not be right, should there be derivs?
		{
			flow.dV(0, 2) = 0.0; 
			flow.dV(1, 2) = 0.0; 
			flow.dV(2, 2) = 0.0; 
			flow.dT(2)    = 0.0; 
			flow.v(2)    = 0.0; 
		}
		break;
		
		default:							// catch boundary conditions yet to be implemented
			cout << "No handler for boundary condition " << ibndcnd << endl;
			cin.get();
		break;
	}
};

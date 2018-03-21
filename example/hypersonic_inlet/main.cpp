#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;
/********************************************************************************
 * Hypersonic inlet grid generator.
 * USAGE
 * Modify
 * Compile (g++ main.cpp)
 * Execute (./a.out)
 * Convert to binary (CONVERT .)
 * Execute (mpiexec -np 8 CFDV hypersonic_inlet)
 *******************************************************************************/


int main()
{
	// works, just add header and section output, test as a whole file generator

	ofstream control("control",ios::out);
	ofstream aero("aero",ios::out);
	ofstream nodes("nodes",ios::out);
	ofstream elements("elements",ios::out);
	ofstream faces("faces",ios::out);

	std::vector<double> Uplate;
	std::vector<double> Uinflow;
	std::vector<double> Udomain;

	std::vector< std::vector< int > > bc;

	Uplate.resize(4);
	Uinflow.resize(4);
	Udomain.resize(4);

	double M      = 20;												// Mach number
	double Tinf   = 1.0;
	double Tplate = (1 + (0.4/2) * pow(M,2));
	double Cv     = 1/1.4/.4/pow(M,2);
	double rhoinf = 1.0;
	double Vinf   = 1.0;

	aero << "M = " << M << endl;
	aero << "Pr = 0.72" << endl;
	aero << "Re = 3000" << endl;
	aero << "Tinf = " << 216.666 << endl;
	aero << "Twall = " << Tplate << endl;
	aero << "Adiabatic = 0" << endl;
	aero << "gamma = 1.4"; 

	Uinflow[0] = rhoinf;
	Uinflow[1] = rhoinf * Vinf;
	Uinflow[2] = 0.0;
	Uinflow[3] = rhoinf * ( Cv * Tinf + 0.5 * pow(Vinf, 2) );

	Uplate[0] = rhoinf;
	Uplate[1] = 0.0;
	Uplate[2] = 0.0;
	Uplate[3] = rhoinf * ( Cv * Tplate );

	Udomain[0] = rhoinf;
	Udomain[1] = rhoinf * Vinf; // was 0.5 Vinf
	Udomain[2] = 0;
	Udomain[3] = rhoinf * ( Cv * Tinf + 0.5 * pow(Vinf,2) );	// was 0.5 Vinf

		// geometry size
	double xnod = 1600;			// 1100
	double ynod = 500;			// 700

	double xele = xnod - 1;
	double yele = ynod - 1;

		// domain size
	double x = 1.6;
	double y = 0.25;

	control << "tmax = 100" << endl;
	control << "itermax = 100000" << endl;
	control << "cfl = 1.0" << endl;
	control << "gmresrestart = 1" << endl;
	control << "gmresiter = 100" << endl;
	control << "ndim = 2" << endl;
	control << "nnod = 4" << endl;
	control << "neqn = 4" << endl;
	control << "nbnod = 2" << endl;
	control << "nface = 4" << endl;
	control << "nodes = " << xnod * ynod << endl;
	control << "elements = " << xele * yele << endl;
	control << "iter = 0" << endl;
	control << "adap = 0";

	double e = 1;
	for (double i=0; i<xnod; i++)
	{
		for (double j=0; j<ynod; j++)
		{
			nodes << e++ << " ";								// node count
			nodes << (i/(xnod-1))*x - 0.1 << " ";			// x coordinate
			nodes << (j/(ynod-1))*y << " ";					// y coordinate
//			nodes << sin((j/(ynod-1))*3.14159-3.14159/2.)*y/2+y << " ";	// y coordinate

			if ( j == 0 && ((i/(xnod-1))*x-0.1) >= 0.0)
			{
				// wall
				for (int k=0; k<4; k++)
				{
					nodes << Uplate[k] << " ";
				}
				nodes << 1;		// wall bc
			}
			else if (i == 0)
			{
				// inflow
				for (int k=0; k<4; k++)
				{
					nodes << Uinflow[k] << " ";
				}
				nodes << -1;		// inflow bc
			}
			else if ( j == ynod-1)
			{
				if (((i/(xnod-1))*x-0.1) >= 0.05)
				{
						// wall
					for (int k=0; k<4; k++)
					{
						nodes << Uplate[k] << " ";
					}
					nodes << 1;		// wall bc
				} else {
						// domain
					for (int k=0; k<4; k++)
					{
						nodes << Udomain[k] << " ";
					}
					nodes << 7;		// top domain
				}
			}
			else
			{
				// domain
				for (int k=0; k<4; k++)
				{
					if (j < 50 && ((i/(xnod-1))*x-0.1) >= 0.0)
					{
						nodes << (j/49)*Udomain[k]+((49-j)/49)*Uplate[k] << " ";
					}
					else if (j > ynod-50 && ((i/(xnod-1))*x-0.1) >= 0.05)
					{
						double jj = ynod - j;
						nodes << (jj/49)*Udomain[k]+((49-jj)/49)*Uplate[k] << " ";
					}
					else
					{
						nodes << Udomain[k] << " ";
					}
				}

				// test for alternate bc's
				if (j == 0)			// wall is caught in first if statement
				{
					nodes << 2;		// reflection
				}
		//		else if ( j == ynod-1)
		//		{
		//			nodes << 7;		// farfield (7, is 1 for test)
		//		}
				else if ( i == xnod-1)
				{
					nodes << 6;		// outflow
				}
				else
				{
					nodes << 0;		// interior domain
				}
			}
			nodes << endl;
		}
	}
	nodes << endl;

	e = 1;

	// Connectivity
	for (int j=1; j<xnod; j++)
	{
		for (int i=1; i<ynod; i++)			// count elements in ynod direction
		{
			elements << e << "\t" << i+(j-1)*ynod << "\t" << i+j*ynod << "\t" << (i+1)+j*ynod << "\t" << (i+1)+(j-1)*ynod << endl;
			e++;
		}
	}
	elements << endl;

	// Attempt at face generation
	
	// inflow
	for (int i=yele-1; i>=0; i--)
	{
		std::vector < int > boundary;
		boundary.push_back( i);					// element number
		boundary.push_back(3);					// face
		boundary.push_back(-1);					// bc
		boundary.push_back(3);					// node0
		boundary.push_back(0);					// node1
		bc.push_back( boundary);
	}
	for (int i=0; i<xele; i++)
	{
		std::vector < int > boundary;
		boundary.push_back( i*yele);			// element number
		boundary.push_back(0);					// face
		if ( (i/(xnod-1))*x - 0.1 < 0)
		{
			boundary.push_back(2);					// bc	
		}
		else
		{
			boundary.push_back(1);					// bc
		}
		boundary.push_back(0);					// node0
		boundary.push_back(1);					// node1
		bc.push_back( boundary);
	}
	for (int i=0; i<yele; i++)
	{
		std::vector < int > boundary;
		boundary.push_back( i+yele*(xele-1));					// element number
		boundary.push_back(1);					// face
		boundary.push_back(6);					// bc
		boundary.push_back(1);					// node0
		boundary.push_back(2);					// node1
		bc.push_back( boundary);
	}
	for(int i=xele-1; i>=0; i--)
	{
		std::vector < int > boundary;
		boundary.push_back( (i+1)*yele-1);					// element number
		boundary.push_back(2);					// face
		if (((i/(xnod-1))*x-0.1) >= 0.05)
			boundary.push_back(1);				// bc			( 7 = farfield, 1 = plate for inlet)
		else
			boundary.push_back(7);
		boundary.push_back(2);					// node0
		boundary.push_back(3);					// node1
		bc.push_back( boundary);
	}
/*
	// need to add corner BC's
	std::vector < int > boundary;
	boundary.push_back(0);					// element number
	boundary.push_back(0);					// face
	boundary.push_back(101);				// bc			( 7 = farfield, 1 = plate for inlet)
	boundary.push_back(0);					// node0
	boundary.push_back(1);					// node1
	bc.push_back( boundary);

	boundary[0] = yele-1;
	boundary[1] = 2;
	boundary[2] = 100;
	boundary[3] = 2;
	boundary[4] = 3;
	bc.push_back( boundary);

	boundary[0] = (xele-1)*yele;
	boundary[1] = 1;
	boundary[2] = 111;
	boundary[3] = 1;
	boundary[4] = 2;
	bc.push_back( boundary);
*/
	for (int i=0; i < bc.size(); i++)
	{
		for (int j=0; j<bc[i].size(); j++)
		{
			faces << bc[i][j] << "\t";
		}
		faces << endl;
	}
	faces << endl;

	control.close();
	aero.close();
	nodes.close();
	elements.close();
	faces.close();

	return 0;
}

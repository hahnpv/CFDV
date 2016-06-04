#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

#include "FDV/Element.h"
#include "FDV/Node.h"

// Output a Tecplot FE data file

struct TecplotOut
{
	TecplotOut(std::vector<Element *> &elements, std::vector<Node *> &nodes, /*int iter,*/ int ndim, /*double time,*/ bool ele_data, int nnod, int neqn)
	{
			// write file
		//	std::ofstream tout(("tecplot_/*to_string<double>(iter)+*/.plt").c_str(),ios::out);
			std::ofstream tout(("tecplot.plt"),ios::out);

				// Write Tecplot header
			tout << "title = \"Output at time = " << time <<" \"" << endl;
			tout << "variables = \"x\", \"y\",";
			if (ndim == 3)
			{
				tout << " \"z\",";
			}
			tout <<	"\"rho\", \"u\", \"v\",";
			if (ndim == 3)
			{
				tout << " \"w\",";
			}		
			tout << "\"T\", \"p\","; 
			if (ele_data)
			{
				tout << "\"s1\", \"s2\", \"s3\", \"s4\", ";
			}
			tout << endl;
			tout << "zone E=" << elements.size() << ", N="<< nodes.size() <<", F=FEBLOCK";
			if (ndim == 3)
			{
				tout << " ET=BRICK";
			}
			tout << endl;
			if (ele_data)
			{
				if ( ndim == 2)
				{
					tout << "VARLOCATION=([8,9,10,11]=CellCentered)" << endl;
				}
				else if (ndim == 3)
				{
					tout << "VARLOCATION=([10,11,12,13]=CellCentered)" << endl;
				}
			}
			tout << endl;
			
				// Set scientific notation w/7 digits of precision
			tout << scientific << endl;
			tout << setprecision(7) << endl;


			// x, y
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->x(0) << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;
			
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->x(1) << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;


			if (ndim == 3)
			{
				for (unsigned int n=0; n<nodes.size(); n++)
				{
					tout << nodes[n]->x(2) << " ";

					if (n%100 == 0 && n != 0) tout << endl;
				}
				tout << endl;
			}

			// rho
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->rho << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

					// u
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->v(0) << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

						// v
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->v(1) << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

			if(ndim == 3)
			{
				for (unsigned int n=0; n<nodes.size(); n++)
				{
					tout << nodes[n]->v(2) << " ";

					if (n%100 == 0 && n != 0) tout << endl;
				}
				tout << endl;
			}

			// temperature
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->T << " ";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;

	// pressure needs to be backed out ... 
			// pressure
			for (unsigned int n=0; n<nodes.size(); n++)
			{
				tout << nodes[n]->p << "\t";

				if (n%100 == 0 && n != 0) tout << endl;
			}
			tout << endl;


			// correlation
			for (unsigned int e=0; e < elements.size(); e++)
			{
				for(unsigned int n=0; n < nnod; n++)
				{
					tout << elements[e]->node[n]->number + 1 << " ";
				}
				tout << endl;
			}

			tout.close();
		}
};
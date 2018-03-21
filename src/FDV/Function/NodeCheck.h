#pragma once

#include <iostream>
using namespace std;

	/// Checks each node for negative or unreal values
template<class T> struct NodeCheck : public unary_function<T, void>
{
	NodeCheck()
	{
	}
	void operator() (T n) 
	{
		if ( n->rho < 0)
		{
			cout << "negative rho in node " << n->number << ", " << n->rho << " ";
			
			for (unsigned int i=0; i < n->x.imax(); i++)
			{
				cout << n->x(i) << " ";
			}
			
			cout << endl;
		}

		if ( n->T < 0)
		{
			cout << "negative T in node " << n->number << ", " << n->T << " "; 
			
			for (unsigned int i=0; i < n->x.imax(); i++)
			{
				cout << n->x(i) << " ";
			}
			
			cout << endl;

//			n->T = 0.;
		}

		if ( n->p < 0)
		{
			cout << "negative P in node " << n->number <<  ", " << n->p<< " "; 
					
			for (unsigned int i=0; i < n->x.imax(); i++)
			{
				cout << n->x(i) << " ";
			}
			
			cout << endl;

//			n->p = 0.;
		}
	}
};


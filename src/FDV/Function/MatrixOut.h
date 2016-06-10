#pragma once

#include <iostream>
using namespace std;

#include "../Element.h"
	/// Dump the matrices in the element
template<class T> struct MatrixOut : public unary_function<T, void>
{
	MatrixOut(std::string fname)
	{
		ofstream mout(fname.c_str(),ios::out);

		for (int e = 0; e < elements.size(); e++)
		{
			mout << "element " << e << endl;
			for (int i = 0; i < 40; i++)
			{
				for (int j = 0; j < 40; j++)
				{
					mout << elements[e]->R(i, j) << ", ";	
				}
				mout << ", " << elements[e]->rhs(i);
				mout << endl;
			}
			mout << endl;
		}

		mout.close();
	}
};


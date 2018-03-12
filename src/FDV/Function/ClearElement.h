#pragma once

#include <iostream>
using namespace std;

	/// Clears the transient data in each element
template<class T> struct ClearElement : public unary_function<T, void>
{
	ClearElement()
	{
	}
	void operator() (T e) 
	{
		e->R.clear();
		e->rhs.clear();
	}
};


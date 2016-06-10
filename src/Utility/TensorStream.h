#pragma once 

#include <ostream>
#include "split.h"
#include "tostring.h"

	///////////////////////////////////////////////////////////////////////////
	/// TENSOR STREAM PRINTING 												///
	///////////////////////////////////////////////////////////////////////////

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 1> & t ) 
{
	for (int i=0; i < t.imax(); i++)
	    os << t(i) << ", ";

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 2> & t ) 
{

	for (int i = 0; i < t.imax(); i++)
	{
		for (int j = 0; j < t.jmax(); j++)
		{
			os << t(i, j) << ", ";
		}
		os << std::endl;
	}

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 3> & t ) 
{

	for (int i = 0; i < t.imax(); i++)
	{
		os << "i = " << i << std::endl;

		for (int j = 0; j < t.jmax(); j++)
		{
			for (int k = 0; k < t.kmax(); k++)
			{
				os << t(i, j, k) << ", ";
			}
			os << std::endl;
		}
		os << std::endl;
	}

    return os;
}

template<class T>
std::ostream& operator << ( std::ostream& os, Tensor<T, 4> & t ) 
{

	for (int i = 0; i < t.imax(); i++)
	{
		for (int j = 0; j < t.jmax(); j++)
		{
			os << "(i, j) = " << i << ", " << j << std::endl;

			for (int k = 0; k < t.kmax(); k++)
			{
				for (int l = 0; l < t.lmax(); l++)
				{
					os << t(i, j, k) << ", ";
				}
				os << std::endl;
			}
			os << std::endl;
		}
	}

	return os;
}



	///////////////////////////////////////////////////////////////////////////
	/// TENSOR SERIALIZATION												///
	///////////////////////////////////////////////////////////////////////////
	// FORMAT:
	// 1. Type
	// 2. ndim
	// 3. imax, jmax, kmax, lmax. Up to ndim
	// 4. data read straight from data, not in any particular indice order
	// kinda gimp since (1) you need T, ndim to instantiate tensor from the start ...
	// so really you only need the data ... but it works. 
template<class T, int ndim>
std::ostream& write ( std::ostream& os, Tensor<T, ndim> & t )
{
	// lame but all I know to do as of right now
	// or encode in Tensor class... can also typeid().name() but can get fugly
	if (typeid(T) == typeid(double))
	{
		os << "double" << " ";
	}
	else if (typeid(T) == typeid(bool))
	{
		os << "bool" << " ";		
	}
	else if (typeid(T) == typeid(int))
	{
		os << "int" << " ";		
	}

	os << ndim;

	if (ndim > 0)
		os << " " << t.imax();
	if (ndim > 1)
		os << " " << t.jmax();
	if (ndim > 2)
		os << " " << t.kmax();
	if (ndim > 3)
		os << " " << t.lmax();

	for (int i = 0; i < t.imax() * t.jmax() * t.kmax() * t.lmax() ; i++)
	{
		os << " " << t.data[i];
	}

	os << endl;

//	cin.get();
	return os;
}

template<class T, int ndim>
std::istream& read ( std::istream& os, Tensor<T, ndim> & t )
{
	std::string line;
	getline(os, line);

	std::vector<std::string> data = split<std::string>(line, " ");

	int count = 0;

	std::string type = data[count++];
/*
	if (type == "double")
	{
		cout << "double" << endl;
	}
	else if (type == "bool")
	{
		cout << "bool" << endl;
	}
	else if (type == "int")
	{
		cout << "int" << endl;
	}
*/
		// should be string_to<int> ?
	int ndimi = string_to<int>(data[count++]);

	int imax = 1, jmax=1, kmax=1, lmax=1;
	if (ndimi > 0)
		imax = string_to<int>(data[count++]);
	if (ndimi > 1)
		jmax = string_to<int>(data[count++]);
	if (ndimi > 2)
		kmax = string_to<int>(data[count++]);
	if (ndimi > 3)
		lmax = string_to<int>(data[count++]);

	for (int i = 0; i < imax * jmax * kmax * lmax ; i++)
	{
		t.ref->data[i] = string_to<T>(data[count++]);
	}

	// FORMAT:
	// 1. Type
	// 2. ndim
	// 3. imax, jmax, kmax, lmax. Up to ndim
	// 4. data read straight from data, not in any particular indice order

	// lame but all I know to do as of right now
	// or encode in Tensor class... can also typeid().name() but can get fugly

//	cin.get();
	return os;
}

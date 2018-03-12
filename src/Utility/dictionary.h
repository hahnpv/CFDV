#pragma once

#include "tostring.h"

// a dictionary 
// might make a singleton ... 
#include <vector>
//template<class T> 
struct dictionary 
{
	dictionary()
	{}

	void add_key(std::vector<std::string>  in)
	{
		if (in.size() >= 2)
		{
			if (in[0].compare("!") != 0 )
			{
				key.push_back(in[0]);
				val.push_back(in[2]);
			}
		}
	}

	void add_key(std::string keyin, std::string valin)
	{
			key.push_back(keyin);
			val.push_back(valin);
	}

	template < class T>
	T get_val(std::string keyval)
	{
		for (unsigned int i=0; i<key.size(); i++)
		{
			if (key[i].compare(keyval) == 0)
			{
				return string_to<T>(val[i]);
			}
		}

		// throw an exception, in the future
		T t(0);
		return t;
	}

	std::vector<std::string> key;
	std::vector<std::string> val;
};

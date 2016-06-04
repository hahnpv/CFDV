#pragma once

#include <sstream>

template <class T>
inline std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

template <class TClass>
TClass string_to(std::string s)
{
	TClass result;
	std::istringstream is( s );
	is >> result;
	return result;
}


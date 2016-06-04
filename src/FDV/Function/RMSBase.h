#pragma once

// RMSErr base
template<class T> struct RMSBase
{
	RMSBase()
	{}

	typedef void result_type;					// result_type for tr1 wrapper.
	virtual result_type operator() (T n) 
	{
    }
};


#pragma once

template<typename T> class Singleton
{
  public:
    static T& Instance()
    {
        static T theSingleInstance;  // assumes T has a protected default constructor
        return theSingleInstance;
    }
};

// example:

/*

class ThermoInstance : public Singleton<ThermoInstance>
{
    friend class Singleton<ThermoInstance>;
    int example_data;
  public:
    int Getexample_data() const {return example_data;}
  protected: 
    ThermoInstance(): example_data(42) {}   // default constructor 
};
 
#define Thermo ThermoInstance::Instance()
 
// This test case should print "42". 
#include <iostream>
int main()
{
    std::cout<< Thermo.Getexample_data()<<std::endl;
	std::cin.get();
    return 0;
}

*/


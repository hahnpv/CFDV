#ifdef __GNUC__
	#include <tr1/functional>			/// ref wrapper (STL TR1 gcc)
	using namespace std::tr1;
#else
	#include <functional>					/// ref wrapper (STL TR1 msvc)
// Floating Point Break on NaN (windows debugging)
//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
#endif

#include "Utility/Tensor.h"

#include <iostream>

using namespace std;

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "MPIBinaryConfiguration.h"

#include "FDV/Function/ClearElement.h"
#include "FDV/Function/NodeUnpack.h"

// Extract data from the elements, faces and specified nodes file
// Start with Tecplot, add Ensight, point extraction, etc.
// Eventually would be nice to do turbulence statistics, etc.

po::variables_map add_program_options(int argc, char * argv[])
{
	// program options
	po::options_description desc("Command Line Parameters:");
	desc.add_options()
	("help", "This help message")
	("path",    po::value<std::string>(), "Path to CFD case")
	("iter",    po::value<int>()->default_value(0), "Iteration")
	("time",    po::value<double>()->default_value(0.0), "Time")
	("output",  po::value<std::string>(), "Output type (Tecplot/Ensight/?)");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || argc <= 1)
	{
		cout << desc << endl;
	}

	return vm;
}

int main(int argc, char *argv[])
{
	std::vector<std::string> args;

	std::vector< Node *> nodes;
	std::vector< Element *> elements;

	po::variables_map vm = add_program_options(argc, argv);

	std::string path;
	std::string nodefile;

	if (vm.count("help") || argc <= 1)
	{
		return -1;
	}

	if (vm.count("path"))
	{
		path = vm["path"].as<std::string>();
		cout << "CFD Path: " << path << endl;

		args.push_back("Extract");
		args.push_back(path);
	}
	else
	{
		cout << "Path to CFD Case not specified!" << endl;
		return -1;
	}

	MPIBinaryConfiguration * config = new MPIBinaryConfiguration(argc, argv, vm, elements, nodes, vm["path"].as<std::string>());

	int iter = vm["iter"].as<int>();
	double t = vm["time"].as<double>();

	// Tecplot
	config->Output(elements, nodes, iter, t, true);

	// Ensight TODO

	return 0;
}

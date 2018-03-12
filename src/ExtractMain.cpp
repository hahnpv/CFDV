#include "Utility/Tensor.h"

#include <iostream>

using namespace std;

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "TecplotOut.h"
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"

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
	("node",    po::value<std::string>(), "Node file")
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

	if (vm.count("node"))
	{
		nodefile = vm["node"].as<std::string>();
		cout << "Node file: " << nodefile << endl;
	}
	else
	{
		cout << "Node file not specified!" << endl;
		return -1;
	}

	LoadBinaryData init(vm, elements, nodes);					/// new binary method
//	init.read_nodes(0, init.key.get_val<int>("nodes") );				// add new node path here.
//	init.read_elements(0, init.key.get_val<int>("elements") );

	int ndim = init.cm["ndim"].as<int>();
	int nnod = init.cm["nnod"].as<int>();
	int neqn = init.cm["neqn"].as<int>();

	for_each(elements.begin(), elements.end(), ClearElement<Element *>());					/// Clear matrices in each element
	for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());								/// extract nodes

	// Tecplot - only need nodes, elements.
	TecplotOut(elements, nodes, /*iter,*/ ndim, /*time,*/ false, nnod, neqn);

	// Ensight - test.
}

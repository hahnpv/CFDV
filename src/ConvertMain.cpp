#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "FDV/Element.h"
#include "FDV/Node.h"

#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
#include "Utility/IO/ASCIItoBinary/SaveBinary.h"
#include "Utility/IO/NewLoad.h"

#include "FDV/Function/NodeUnpack.h"

#include "Utility/MPI_Breakdown.h"

using namespace std;

int main(int argc, char *argv[])
{
	std::vector<std::string> args;
	for (int i=0; i<argc; i++)
	{
		args.push_back( string(argv[i]) );
	}

	std::vector< Node *> nodes;
	std::vector< Element *> elements;

	cout << "reading input ... " << endl;
//	LoadBinaryData init(args, elements, nodes);								/// load using new input method
	NewLoad init(args,elements,nodes);
	int nnod = init.key.get_val<int>("nnod");							/// Number of nodes per finite element
	int neqn = init.key.get_val<int>("neqn");							/// Number of equations per node
	int ndim = init.key.get_val<int>("ndim");
	std::string path = init.key.get_val<std::string>("path");

	cout << "breaking down MPI ... " << endl;
	MPI_Breakdown::clear(path);
	for (int i = 1; i <= 8; i++)
	{
		MPI_Breakdown breakdown( elements, nodes, neqn, nnod, path, i);
	}

	for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());

	cout << "saving binary data ... " << endl;
	if (ndim == 3)
	{
		SaveBinaryData<ele_3d_t, node_3d_t> save(elements, nodes, ndim, neqn, nnod);
		save.write_elements(path + "//binaryelements");
		save.write_nodes(path + "//binarynodes");
	}
	else if (ndim == 2)
	{
		SaveBinaryData<ele_2d_t, node_2d_t> save(elements, nodes,  ndim, neqn, nnod);
		save.write_elements(path + "//binaryelements");
		save.write_nodes(path + "//binarynodes");
	}
	else if (ndim == 1)
	{
		SaveBinaryData<ele_1d_t, node_1d_t> save(elements, nodes,  ndim, neqn, nnod);
		save.write_elements(path + "//binaryelements");
		save.write_nodes(path + "//binarynodes");
	}
	cout << "done!" << endl;
//	cin.get();
	return 0;
};

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "FDV/Element.h"
#include "FDV/Node.h"

#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
#include "Utility/IO/ASCIItoBinary/SaveBinary.h"
#include "FDV/Function/NodeUnpack.h"
#include "Utility/IO/AsciiLoad.h"

#include "Utility/MPI_Breakdown.h"

using namespace std;

/// CONVERT exists to convert text node/element correlations to binary
int main(int argc, char *argv[])
{
	std::vector<std::string> args;
	for (int i=0; i<argc; i++)
	{
		args.push_back( string(argv[i]) );
		cout << "argument " << string(argv[i]) << endl;
	}

	std::vector< Node *> nodes;
	std::vector< Element *> elements;

	cout << "reading input ... " << endl;

	AsciiLoad init(args,elements,nodes);
	int nnod  = (int)(init.key.find("nnod")->second);
	int neqn  = (int)(init.key.find("neqn")->second);
	int ndim  = (int)(init.key.find("ndim")->second);
	std::string path = init.path;

	for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());
	cout << "saving binary data ... " << endl;
	std::string savepath = path;
	boost::filesystem::create_directory(savepath);

	// SaveBinaryData makes implicit MPI assumptions. May revisit
	int rc = MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// this is ripped from MPIBinaryConfiguration. Might find a way to make this generic/accessible
	SaveBinaryData * save = new SaveBinaryData(elements, nodes, ndim, neqn, nnod, rank, savepath);
	save->write_nodes(0,       nodes.size());
	save->write_elements(0, elements.size());
	save->write_adap();

	MPI_Breakdown::clear(savepath);
	for(int i=1; i<17;i++)
		MPI_Breakdown(elements, nodes, neqn, nnod, savepath, i);

	cout << "done!" << endl;

	MPI_Finalize();

	return 0;
};

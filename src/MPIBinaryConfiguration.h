#pragma once

#include "Utility/IO/MPI_TecplotOut.h"
#include "FDV/Function/MPI_RMSError.h"
#include "Utility/Solvers/MPI_GMRES.h"
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
#include "Utility/IO/ASCIItoBinary/SaveBinary.h"
#include "Utility/MPI_Breakdown.h"

#include <iomanip>
#include "boost/filesystem.hpp"				// includes all needed Boost.Filesystem declarations
#include "boost/filesystem/fstream.hpp"		// filestreams, makes avail to other classes and may cause problems?

#include <functional> // std::ref

#include "boost/program_options.hpp"
namespace po = boost::program_options;
using namespace std;

struct MPIBinaryConfiguration
{
	MPIBinaryConfiguration(int argc, char * argv[], po::variables_map & cmdline, std::vector<Element *> & elements, std::vector<Node *> & nodes, std::string path)
		: path(path)
	{
		///// NEW
		path = cmdline["path"].as<std::string>();

		// input file variable map
		po::options_description controlinput("Configuration File Options");
		controlinput.add_options()
			("tmax", po::value<double>(), "Maximum flowthrough time [s]")
			("itermax", po::value<int>(), "Maximum number of iterations")
			("cfl", po::value<double>(), "Prescribed CFL number")
			("gmresrestart", po::value<int>(), "GMRES: number of restarts")
			("gmresiter", po::value<int>(), "GMRES: number of iterations")
			("ndim", po::value<int>(), "Number of nodes")
			("nnod", po::value<int>(), "Number of dimensions")
			("neqn", po::value<int>(), "Number of equations")
			("nbnod", po::value<int>(), "Number of boundary nodes")
			("nface", po::value<int>(), "Number of faces")
			("nodes", po::value<int>(), "Number of nodes")
			("elements", po::value<int>(), "Number of elements")
			("iter", po::value<int>(), "Starting iteration number")
			("adap", po::bool_switch()->default_value(false), "Use adaption file to specify hanging nodes")
			("axi", po::bool_switch()->default_value(false), "Axisymmetric (2D only)")
			;

		std::string fname = path + "//control";
		std::ifstream control_file(fname.c_str(), std::ifstream::in);
		po::store(po::parse_config_file(control_file, controlinput), cm);
		control_file.close();
		po::notify(cm);

		po::options_description aeroinput("Configuration File Options");
		aeroinput.add_options()
			("M", po::value<double>(), "Freestream Mach number")
			("Pr", po::value<double>(), "Prandtl Number")
			("Re", po::value<double>(), "Freestream Reynolds number")
			("Tinf", po::value<double>(), "Freestream Temperature")
			("Twall", po::value<double>(), "Wall Temperature")
			("Adiabatic", po::value<bool>(), "Adiabatic Wall flag")
			("gamma", po::value<double>(), "Ratio of Specific Heats")
			;

		fname = path + "//aero";
		std::ifstream aero_file(fname.c_str(), std::ifstream::in);
		po::store(po::parse_config_file(aero_file, aeroinput), am);
		aero_file.close();
		po::notify(am);

		int rc = MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		// TODO combine control+aero into one file with [headers]

		// set up nnod, neqn, ndim
		nnod = cm["nnod"].as<int>();							/// Number of nodes per finite element
		neqn = cm["neqn"].as<int>();							/// Number of equations per node
		ndim = cm["ndim"].as<int>();							/// Number of dimensions 
		int nbnod = cm["nbnod"].as<int>();
		gmres_iter = cm["gmresiter"].as<int>();					/// Number of GMRES iterations
		gmres_rest = cm["gmresrestart"].as<int>();				/// Number of GMRES restarts

		LoadBinaryData load(elements, nodes, path, nnod, neqn, ndim, nbnod);	// TODO move this into config						/// new binary method

		// read_mpi();
		MPI_Breakdown::read(path, rank, size, nodemin, nodemax, elemin, elemax, nedgenode, nedgenode_left, max, node_iterator_min, node_iterator_max);

		Thermo & thermo = Thermo::Instance();
		thermo.set(am["gamma"].as<double>(),
			am["M"].as<double>(),
			am["Re"].as<double>(),
			am["Pr"].as<double>(),
			am["Tinf"].as<double>(),
			am["Twall"].as<double>(),
			am["Adiabatic"].as<bool>());

		load.read_nodes(nodemin, nodemax);		
		load.read_elements(elemin, elemax);	

		cout << "nodes.size: " << nodes.size() << endl;

		NodeIteratorStart = nodes.begin();
		NodeIteratorEnd = nodes.begin() + (node_iterator_max - node_iterator_min) + 1;

		offset = nodemin;
		nodesize = (nodes.size() - nedgenode) * neqn;

		if (cm["adap"].as<bool>())
		{
			load.read_adap();
		}

		cout << "rank:          " << rank << endl;
		cout << "element range: " << elemin << " " << elemax << " " << elements.size() << endl;
		cout << "node range:    " << nodemin << " " << nodemax << " " << nodes.size() << endl;
		cout << "nedgenode:     " << nedgenode << " left: " << nedgenode_left << endl;
		cout << "node iter:     " << 0 << " to " << node_iterator_max - node_iterator_min + 1 << endl;//" which is " 
		cout << "offset:        " << offset << endl;
		cout << "node size:     " << nodesize << endl;

		gmres = new MPI_GMRES( gmres_rest, gmres_iter, nodesize, neqn, nnod, ndim, size, rank, offset);
	}
	
	~MPIBinaryConfiguration()
	{
		MPI_Finalize();
		delete gmres;
	}

	void reset_gmres()
	{
		delete gmres;
		gmres = new MPI_GMRES( gmres_rest, gmres_iter, nodesize, neqn, nnod, ndim, size, rank, offset);
	}

	void RMSErr(double t, int iter)
	{
		MPI_RMSError<Node *> rmserr("RMSError.txt", t, iter, neqn);
		for_each(NodeIteratorStart, NodeIteratorEnd, ref(rmserr));
	}
	void Output(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter, double time, bool ele_data)
	{
		MPI_TecplotOut(elements, nodes, iter, ndim, time, ele_data , rank, offset, size, nedgenode_left, nnod, neqn);
	}

	void Save(std::vector<Element *> &elements, std::vector<Node *> &nodes, int iter)
	{
		std::string savepath = path + "//" + to_string<int>(iter) + "//";
		boost::filesystem::create_directory(savepath);

		SaveBinaryData * save = new SaveBinaryData(elements, nodes, ndim, neqn, nnod, rank, savepath);
		save->write_nodes(nodemin, nodemax);
		save->write_elements(elemin, elemax);
		save->write_adap();

		MPI_Breakdown::clear(savepath);
		for(int i=1; i<17;i++)
			MPI_Breakdown(elements, nodes, neqn, nnod, savepath, i);
	}

	int size;
	int offset;
	int nedgenode;
	int nedgenode_left;
	int nodesize;
	int nodemin, nodemax;
	int elemin, elemax;
	int node_iterator_min, node_iterator_max; // allows functionalizatoin, don't like this.
	std::string path;

	int max;				// int equiv of max iterator

	int gmres_iter, gmres_rest;

	int rank;

	int nnod;
	int neqn;
	int ndim;
	std::vector<Node *>::iterator NodeIteratorStart;
	std::vector<Node *>::iterator NodeIteratorEnd;
	GMRESbase * gmres;

public:
	po::variables_map am;
	po::variables_map cm;
};

#include "mpi.h"

//#ifdef __GNUC__
//	#include <tr1/functional>			/// ref wrapper (STL TR1 gcc)
//	using namespace std::tr1;
//#else
	#include <functional>					/// ref wrapper (STL TR1 msvc)
// Floating Point Break on NaN (windows debugging)
///unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
//#endif

#include <algorithm>					/// for_each

	// Configuration
#include "MPIBinaryConfiguration.h"

	// Base Classes
#include "FDV/Element.h"
#include "FDV/Node.h"
#include "FDV/FDVGalerkin.h"

	// Dirichlet BC
#include "FDV/Function/ApplyBC.h"
#include "FDV/Function/EnforceBC.h"

	// CFD Utility
#include "FDV/Function/CalcLength.h"
#include "FDV/Function/ClearElement.h"
#include "FDV/Function/NodeUnpack.h"
#include "FDV/Function/NodeCheck.h"
#include "FDV/Function/MPI_CFL.h"
#include "FDV/Function/FDVParam.h"

	// Input
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"

	// Sensors
#include "Utility/Sensor.h"

	// Non-CFD
#include <boost/timer.hpp>
#include "Utility/Timer.h"

#include "Utility/PrecisionTimer.h"

	// Mesh Refinement
#include "adap/ElementAssociation.h"
#include "adap/MeshRefine2D.h"
#include "adap/RefineListAdaptive.h"				// attempt at true adaption based on s1 parameter

	// Test
#include "FDV/Function/AxisymmetricFlow.h"

#include "boost/program_options.hpp"
namespace po = boost::program_options;

using namespace std;

po::variables_map add_program_options(int argc, char * argv[])
{
	// program option variable map
	po::options_description desc("Command Line Parameters:");
	desc.add_options()
		("path", po::value<std::string>(),				  "Path to CFD case")
		("adap", po::bool_switch()->default_value(false), "Grid Adaptation flag")
		;

	po::positional_options_description p;
	p.add("path", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (argc <= 1)
	{
		cout << desc << endl;
		cout << "The Flowfield Dependent Variation (FDV) method is a numerical scheme" << endl
			<< "built to simulate flows characterized by multiple speeds, multiple " << endl
			<< "physical phenomena, and large variations in flow variables." << endl
			<< "Written by Dr. Philip Hahn, derived from Dr. T.J. Chung's textbook " << endl
			<< "'Computational Fluid Dynamics, 2nd Ed.' and courses CFD I,II,III at " << endl
			<< "The Univerisity of Alabama in Huntsville." << endl
			<< endl
			<< "hahnpv@gmail.com" << endl
			<< endl
			<< "'You will never forget.' -Dr. Chung" << endl;
	}

	return vm;
}

	//	3. Find unifications between 2D and 3D if any
	// Still some hardcoded stuff that needs to be generalized.
int main(int argc, char *argv[])
{
	po::variables_map vm = add_program_options(argc, argv);
	if (argc <= 1)
		return 0;

	double t = 0;
	double dt = 0;

	std::vector< Node *> nodes;
	std::vector< Element *> elements;

	MPIBinaryConfiguration * config = new MPIBinaryConfiguration(argc, argv, vm, elements, nodes, vm["path"].as<std::string>());

	int iter   = config->cm["iter"].as<int>();
	double cfl = config->cm["cfl"].as<double>();						/// Prescribed CFL number
	int itermax = config->cm["itermax"].as<int>();						/// Maximum number of time iterations
	int nnod  = config->cm["nnod"].as<int>();							/// Number of nodes per finite element
	int neqn  = config->cm["neqn"].as<int>();							/// Number of equations per node
	int ndim  = config->cm["ndim"].as<int>();							/// Number of dimensions 
	int nface = config->cm["nface"].as<int>();							/// Number of faces
	int nbnod = config->cm["nbnod"].as<int>();							/// Number of boundary nodes

	cout << "nodes/eles: " << nodes.size() << " " << elements.size() << endl;

	CalcLength<Element *> cl(ndim, nnod);
	for_each(elements.begin(), elements.end(), ref( cl));		// calculate characteristic length

	Sensor * sensors = new Sensor(vm["path"].as<std::string>() + "//sensors", ndim, nodes);

	// Select equations to solve
	NavierStokes * NS;
	if (ndim == 1)
		NS = new NavierStokes1D;
	else if (ndim == 2)
		NS = new NavierStokes2D;
	else
		NS = new NavierStokes3D;

	// can move somewhere else? its a mess here
	IterationTimer timer(config->rank);
	ofstream tout("timeloop.csv",ios::out);
	cout << "time loop" << endl;

PrecisionTimer timer_l;
timer_l.start();

	/// MESH REFINEMENT ///
	if (vm["adap"].as<bool>())
	{
		cout << "Performing grid refinement..." << endl;
		for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());								/// extract nodes
		for_each(nodes.begin(), nodes.end(), NodeCheck<Node *>());								/// check for invalid nodal values
		std::vector<int> elelist;			// elements to be refined
		std::vector<int> refine_level;		// breakpoints for each level refine
		AdaptiveRefineList<Element, Node>(elelist, refine_level, elements, nodes, nnod, nface);
		MeshRefine2D<Element, Node>(elements, nodes, elelist, refine_level, ndim, neqn, nnod, iter, nface);
		cout << "Done adapting mesh" << endl;
		for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());								/// extract nodes
		for_each(nodes.begin(), nodes.end(), NodeCheck<Node *>());								/// check for invalid nodal values
		config->Save(elements, nodes, iter);													// TODO add updated config files?
		return 0;
	}

	config->Save(elements, nodes, iter);														// TODO add updated config files?

	/// CFD Integration ///
	for (; iter < itermax; iter++)
	{
		if (config->rank == 0)
			cout << " t = " << t << ", iteration " << iter << endl;

		timer.reset();

		for_each(elements.begin(), elements.end(), ClearElement<Element *>());					/// Clear matrices in each element

		for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());								/// extract nodes

		for_each(nodes.begin(), nodes.end(), NodeCheck<Node *>());								/// check for invalid nodal values

		if (config->rank == 0)	sensors->update(iter);											// FIXME: this means rank 0 has all nodes!?

		if (iter % 1000 == 0) config->Save(elements, nodes, iter);								// TODO add updated config files?

		CFL<Element *>(elements, cfl, dt);														/// Update CFL number

		if (config->rank == 0) cout << "dt: " << dt << endl;

		if (iter % 100 == 0 )	config->Output(elements, nodes, iter, t, true);					/// Tecplot output, deprecate for Save

		FDVGalerkin<Element *> fdv(nnod, neqn, ndim, nbnod, dt, NS);
		for_each(elements.begin(), elements.end(), ref(fdv));

		if (config->cm["axi"].as<bool>())
		{
			cout << "Processing 2D mesh as axisymmetric" << endl;									// until you validate
			for_each(elements.begin(), elements.end(), AxisymmetricFlow<Element *>(neqn, nnod));	/// Axisym 2pi
		}

		for_each(elements.begin(), elements.end(), ApplyBC<Element *>(neqn, nnod, ndim, nbnod));			/// FIXME this doesn't work at all

		config->gmres->iterate(elements);														/// Solve via GMRES
		config->gmres->update(nodes);
		config->RMSErr(t, iter);																/// Calculate residuals FIXME for adap grids

		t += dt;

		timer_l.stop();
		if(config->rank == 0)
			cout << "time: " << timer_l.read() << endl;
		timer.check(iter, t, dt);
	}

	config->Save(elements, nodes, iter);								// TODO add updated config files?

	delete sensors;
	delete config;
	delete NS;

	return 0;
}

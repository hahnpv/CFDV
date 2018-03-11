#include "mpi.h"

#ifdef __GNUC__
	#include <tr1/functional>			/// ref wrapper (STL TR1 gcc)
#else
	#include <functional>				/// ref wrapper (STL TR1 msvc)
#endif

#include <algorithm>					/// for_each

	// Configuration
#include "ConfigurationBase.h"
//#include "MPIConfiguration.h"
#include "MPIBinaryConfiguration.h"

	// Base Classes
#include "FDV/Element.h"
#include "FDV/Node.h"
#include "FDV/FDVGalerkin.h"

	// Dirichlet BC
//#include "FDV/Function/ApplyBC.h"
#include "FDV/Function/EnforceBC.h"

	// CFD Utility
#include "FDV/Function/CalcLength.h"
#include "FDV/Function/ClearElement.h"
#include "FDV/Function/NodeUnpack.h"
#include "FDV/Function/NodeCheck.h"
#include "FDV/Function/MPI_CFL.h"
#include "FDV/Function/FDVParam.h"
// #include "FDV/Function/MatrixOut.h"
// #include "LUFactorization.h"
//#include "FDV/Function/Calc_Coeff.h"

	// Input
#include "Utility/IO/ASCIItoBinary/LoadBinary.h"
//#include "Utility/IO/NewLoad.h"
// #include "LoadFromFile.h"
// #include "LoadFromTecplot.h"

	// Non-CFD
#include <boost/timer.hpp>
#include "Utility/Timer.h"

#include "Utility/PrecisionTimer.h"

	// Mesh Refinement
#include "adap/ElementAssociation.h"
#include "adap/MeshRefine2D.h"
// #include "MeshUnrefine2D.h"
#include "adap/RefineListAdaptive.h"				// attempt at true adaption based on s1 parameter
// #include "RefineListProgrammed.h"			// programmed test of element refinement

	// Test
#include "FDV/Function/AxisymmetricFlow.h"

	// Floating Point Break on NaN (windows debugging)
//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);


using namespace std;

	//	3. Find unifications between 2D and 3D if any
	// Still some hardcoded stuff that needs to be generalized.
	
int main(int argc, char *argv[])
{
	// deprecate args and replace with program_options, since order could become convoluted.
	std::vector<std::string> args;
	for (int i=0; i<argc; i++)
	{
		args.push_back( string(argv[i]) );
	}

	std::vector< Node *> nodes;
	std::vector< Element *> elements;

	LoadBinaryData init(args, elements, nodes);							/// new binary method

	double t    = 0;
	double dt   = 0;
	int    iter = init.key.get_val<int>("iter");
	double cfl  = init.key.get_val<double>("cfl");						/// Prescribed CFL number
	int itermax = init.key.get_val<int>("itermax");						/// Maximum number of time iterations
	int nnod = init.key.get_val<int>("nnod");							/// Number of nodes per finite element
	int neqn = init.key.get_val<int>("neqn");							/// Number of equations per node
	int ndim = init.key.get_val<int>("ndim");							/// Number of dimensions 
	int nface = init.key.get_val<int>("nface");							/// Number of faces
	int nbnod = init.key.get_val<int>("nbnod");							/// Number of boundary nodes
//	ele_t eletype = init.key.get_val<ele_t>("ele_t");					/// element type

	ConfigurationBase *	config = new MPIBinaryConfiguration(argc, argv, elements, nodes, init.key, init);	

	cout << "nodes/eles: " << nodes.size() << " " << elements.size() << endl;

	cout << "calc length" << endl;
	CalcLength<Element *> cl(ndim, nnod);
	for_each(elements.begin(), elements.end(), std::tr1::ref( cl));		// calculate characteristic length

	// try to move FDV up here, N-S selector struct?

	// can move somewhere else? its a mess here
	IterationTimer timer(config->rank);

	ofstream tout("timeloop.csv",ios::out);

	cout << "time loop" << endl;


PrecisionTimer timer_l;
timer_l.start();

		// Time Loop
	for (; iter < itermax; iter++)
	{
		cout << " t = " << t << ", iteration " << iter << endl;

		timer.reset();

		for_each(elements.begin(), elements.end(), ClearElement<Element *>());					/// Clear matrices in each element
		
		if(iter==1)
		{ 
			/// MESH REFINEMENT ///
			std::vector<int> elelist;			// elements to be refined
			std::vector<int> refine_level;		// breakpoints for each level refine
			AdaptiveRefineList<Element, Node>(elelist, refine_level, elements, nodes, nnod, nface);
			MeshRefine2D<Element, Node>(elements, nodes, elelist, refine_level, ndim, neqn, nnod, iter, nface);
			cout << "Done adapting mesh" << endl;
			/// MESH REFINEMENT ///
		}
		
		for_each(nodes.begin(), nodes.end(), NodeUnpack<Node *>());								/// extract nodes

		for_each(nodes.begin(), nodes.end(), NodeCheck<Node *>());								/// check for invalid nodal values
		
		if (iter == 1)
		{
			config->Save(elements, nodes, iter, init.key);						// Only saves nodes
			cout << "nele: " << elements.size() << endl;
			cout << "nnod: " << nodes.size() << endl;
		}
		
//		if (iter % 10 == 0)	config->Save(elements, nodes, iter, init.key);						// Only saves nodes

		CFL<Element *>(elements, cfl, dt);														/// Update CFL number

		if (config->rank == 0) cout << "dt: " << dt << endl;

		//if (iter % 25 == 0 ) 
		config->Output(elements, nodes, iter, t, true);					/// Tecplot output, deprecate for Save

		/* 
		// extrap consvar within elements for "6" BC
		// 2D quad only FIXME does not work for tri
//		cout << "extrap" << endl;
		for (int e = 0; e < elements.size(); e++)
		{
//			cout << "e = " << e << endl;
			for (int n = 0; n < elements[e]->face.size(); n++)
			{
				if (elements[e]->face[n]->bc == 6)
				{
				//	cout << "found '6' BC in " << e << endl;
					int n0 = elements[e]->face[n]->n(0);
					int n1 = elements[e]->face[n]->n(1);

					int bc0 = n0 - 1;
					int bc1 = n1 + 1;

					if (bc0 < 0)
						bc0 = 3;

					if (bc1 > 3)
						bc0 = 0;


				//	cout << "nodes / bcs: " << n0 << " " << n1 << " " << bc0 << " " << bc1 << endl;
					for (int j = 0; j < neqn; j++)
					{
						elements[e]->node[ n0 ]->U( j) =  elements[e]->node[ bc0 ]->U( j);
						elements[e]->node[ n1 ]->U( j) =  elements[e]->node[ bc1 ]->U( j);
					}
				//	cout << "done swapping consvar" << endl;
				}
			}
		}
		*/

		if (ndim == 3)
		{
			FDVGalerkin<Element *, NavierStokes3D> fdv( nnod, neqn, ndim, nbnod, dt);
			for_each(elements.begin(), elements.end(), std::tr1::ref(fdv));
		}
		else if (ndim == 2)
		{
			FDVGalerkin<Element *, NavierStokes2D> fdv( nnod, neqn, ndim, nbnod, dt);
			for_each(elements.begin(), elements.end(), std::tr1::ref(fdv));
		}
		else if (ndim == 1)
		{
			FDVGalerkin<Element *, NavierStokes1D> fdv( nnod, neqn, ndim, nbnod, dt);
			for_each(elements.begin(), elements.end(), std::tr1::ref(fdv));
		}

//		for_each(elements.begin(), elements.end(), AxisymmetricFlow<Element *>(neqn, nnod));	/// Axisym 2pi

		config->gmres->iterate(elements);														/// Solve via GMRES
		config->gmres->update(nodes);

		for_each(nodes.begin(), nodes.end(), EnforceBC<Node *>(ndim));							/// Enforce Dirichlet boundary conditions
//		config->RMSErr(t, iter);																/// Calculate residuals FIXME for adap grids
		
		t += dt;
		
		timer_l.stop();
		cout << "time: " << timer_l.read() << endl;
		timer.check(iter, t, dt);
	}

	delete config;
    
	MPI_Finalize(); 
	
	return 0;
}

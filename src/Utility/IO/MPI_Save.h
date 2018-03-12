#pragma once

#include <iomanip>
#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations
#include "boost/filesystem/fstream.hpp"   // filestreams, makes avail to other classes and may cause problems?

#include "Thermo.h"

	// Works insofar as it successfully restores a grid
	// Residuals are close but DO NOT MATCH.
	// Probably due to sig figs or something.
	// mess with later see if there is a double output or at least scientific to 8+ places.
struct MPI_Save
{
	MPI_Save( std::vector< Element * > & elements, 
		  std::vector< Node *> & nodes,
		  int ndim,
		  int nnod, 
		  int iter,
		  int rank, 
		  dictionary & key
		)
	{
		Thermo & thermo = Thermo::Instance();

		std::string myRank = to_string<int>(rank);

		// 1. create a save directory
		// 2. copy control / nondim info
		// 3. generate new nodes/elements/faces/adap

		// NOTE: will need to tweak control, not just copy, if we go from non-adap to adap.
		// Also Note: as we convert to new methodology, remove inconsistencies in numbering

//		boost::filesystem::path orig_path( "input//3dtestlagrange" );			// temp eventually dump working directory in db
		boost::filesystem::path orig_path( "input//3dtestlagrange" );			// temp eventually dump working directory in db
		boost::filesystem::path new_path( orig_path / ("restore_pt" + to_string<int>(iter)) );	
		boost::filesystem::create_directory( new_path );

		boost::filesystem::copy_file( orig_path / "aero", new_path / "aero" );				// copy aero file
//		boost::filesystem::copy_file( orig_path / "control", new_path / "control" );		// copy control file
		
		boost::filesystem::ofstream contfile( new_path / "control" / myRank);
		boost::filesystem::ofstream nodefile( new_path / "nodes" / myRank);
		boost::filesystem::ofstream elefile(new_path / "elements" / myRank);
		boost::filesystem::ofstream facefile(new_path / "faces" / myRank);
		boost::filesystem::ofstream adapfile(new_path / "adap" / myRank);

		boost::filesystem::ofstream rankfile(new_path / "rank" / myRank);

		// build nodes
		// essentially replicating NodeUnpack, if we verify we always call after NodeUnpack
		// we could just snag rho, v, T...
		nodefile << setprecision(9) << scientific;
		for (int i = 0; i < nodes.size(); i++)
		{
			Tensor<double, 1> consvar(ndim+2);	// seting == U corrupts u? wtf? problem w/ operator=()?
			consvar(0) = nodes[i]->U(0);

			double vsqr = 0;
			for (unsigned int j=0; j < ndim; j++)
			{
				consvar(j+1) = nodes[i]->U(1+j) / nodes[i]->U(0);
				vsqr += (consvar(j+1) * consvar(j+1));
			}

			consvar(ndim+1) = nodes[i]->U(ndim+1) / nodes[i]->U(0);
			consvar(ndim+1) -= (0.5 * vsqr);
			consvar(ndim+1) /= thermo.Cv;

			nodefile << i+1 << " "
					 << nodes[i]->x << " "
					 << consvar << endl;
		}

		// build elements
		for (int i = 0; i < elements.size(); i++)
		{
			elefile << i+1 << " ";
			for (int j = 0; j < nnod; j++)
			{
				elefile << elements[i]->node[j]->number + 1 << " ";
			}
			elefile << endl;
		}

		// build faces
		for (int i = 0; i < elements.size(); i++)
		{
			if ( elements[i]->face.size() != 0)
			{
				for (int j = 0; j < elements[i]->face.size(); j++)
				{
					facefile << i << " 0 " << elements[i]->face[j]->bc << " " << elements[i]->face[j]->n << endl;
				}
			}
		}

		// build adap
		for (int i = 0; i < elements.size(); i++)
		{
			if ( elements[i]->adap || elements[i]->refine_level > 0)
			{
				// adap flag is redundant unless we wind up doing all eles and output false/identity info.
				adapfile << i << " " << elements[i]->adap << " " << elements[i]->refine_level << " ";
		//		write( adapfile, elements[i]->Hnm);
				cout << "Write is not avialable yet" << endl;
			}
		}	

		// build control file
		contfile << "! control file" << endl;
		contfile << "itermax = " << key.get_val<int>("itermax") << endl;
		contfile << "cfl = " << key.get_val<int>("cfl") << endl;
		contfile << "gmresrestart = " << key.get_val<int>("gmresrestart") << endl;
		contfile << "gmresiter = "	  << key.get_val<int>("gmresiter")    << endl;
		contfile << "ndim = " << key.get_val<int>("ndim") << endl;
		contfile << "nnod = " << key.get_val<int>("nnod") << endl;
		contfile << "neqn = " << key.get_val<int>("neqn") << endl;
		contfile << "adap = " << key.get_val<int>("adap") << endl;
		contfile << "nodes = " << key.get_val<int>("nodes") << endl;
		contfile << "elements = " << key.get_val<int>("elements") << endl;
		contfile << "iter = " << iter << endl;

		// build rank file (specifics like node chopping, etc.)

		contfile.close();
		nodefile.close();
		elefile.close();
		facefile.close();
		adapfile.close();
	}

};
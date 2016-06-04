#pragma once

// new universal loader
// file		- contents
// aero		- as-is
// control	- as-is
// nodes	- nodal data
// elements	- element data
// faces	- face data
// adap		- adaptive meshing data, if needed, put adap flag in control.

#include "Load.h"
#include "FDV/Thermo.h"
#include "FDV/Element.h"
#include "Utility/split.h"
#include "Utility/tostring.h"
#include "Utility/dictionary.h"

	/// Load from file a CFD solution
	/// Restoring from save WORKS except when restoring after first refine going into second refine
	/// not sure why. Might be a precision issue?
	/// but if not refining residuals are identical for 4+ iterations
	/// preserves appropriate adap structure ... so must be precision?
struct NewLoad : public Load
{
	NewLoad( std::vector<std::string> & args,
			 std::vector< Element * > & elements, 
			 std::vector< Node *> & nodes 
			 )
			 : Load(args, elements, nodes)
	{
	
		// args points to DIRECTORY instead of FILES.
		// Load will create bad paths dont use them

		std::string path;
		if (args.size() < 2)			// catch debug instantiation
		{
			path = "..//input//backstep_refined_new";
		//	path = "input//3dtest_new";
	//		path = "..//input//smallbackstep";
		}
		else
		{
			path = args[1];
		}

		key.add_key("path", path);

		cout << "path: " << path << endl;
		cout << "aeropath: " << path + "//aero" << endl;

		aero(path + "//aero");
		control(path + "//control");

		// set up nnod, neqn, ndim
		nnod = key.get_val<int>("nnod");
		neqn = key.get_val<int>("neqn");
		ndim = key.get_val<int>("ndim");
		nbnod = key.get_val<int>("nbnod");

		cout << "nbnod: " << nbnod << endl;

		// set up node/element size here eventually
		cout << "resiszing nodes" << endl;	
		nodes.resize( key.get_val<int>("nodes") );
		cout << "resizing elements" << endl;
		elements.resize( key.get_val<int>("elements") );

		cout << "adding thremo" << endl;
		// need to set up thermo before nodes
		add_thermo();	

		// Load nodes, elements, faces. Add parameters for min and max. 
		// Call MPIConfiguration prior to this and use parameters for our rank.
		// i
		cout << "loading nodes" << endl;
		loadNodes(path + "//nodes");
		cout << "loading elements" << endl;
		loadElements(path + "//elements");
		cout << "loading faces" << endl;
		faces(path + "//faces");

		if ( key.get_val<bool>("adap") )
		{
			cout << "adap" << endl;
			adap(path + "//adap");
		}
	}
	
	void aero(std::string path)
	{
		read_file(path);
	};

	void control(std::string path)
	{
		read_file(path);
	};

	void loadNodes(std::string path)
	{
		std::string input;
		ifstream file(path.c_str());
		while (!file.eof())
		{
			getline(file,input);
			add_node( split<double>(input," \t") );
		}
	};
	void loadElements(std::string path)
	{
		std::string input;
		ifstream file(path.c_str());
		while (!file.eof())
		{
			getline(file,input);
			add_element( split<int>(input," \t") );
		}
	};
	void faces(std::string path)
	{
		std::string input;
		ifstream file(path.c_str());
		while (!file.eof())
		{
			getline(file,input);
			add_face( split<int>(input," \t") );
		}
	};
	void adap(std::string path)
	{
		std::string input;
		ifstream file(path.c_str());
		while (!file.eof())
		{
			getline(file,input);
//			add_adap( split<double>(input," \t") );
			std::vector<double> data = split<double>(input, " \t");
			add_adap( data);
		}
	};

	void add_adap(std::vector<double> & data)
	{
		// 0 - element
		// 1 - adap flag
		// 2 - tensor data
		if (data.size() == 0)
			return;

		int elno = (int)data[0];
		bool adap = (bool)data[1];
		int refine_level = (int)data[2];

		elements[ elno ]->adap = adap;
		elements[ elno ]->refine_level = refine_level;

		std::string datastring = "double ";
		for (int i = 4; i < data.size(); i++)			// this all is hackish, get working then fix.
		{												// might make more sense to make a string ver. of read
			datastring += to_string<double>(data[i]);	// and have the stream version call it once it stringsplits
			datastring += " ";
		}
		std::istringstream is(datastring);

		read( is, elements[ elno ]->Hnm);
	}
};

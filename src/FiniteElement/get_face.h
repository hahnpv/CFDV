#pragma once

#include "../Utility/Tensor.h"
#include "../FDV/Element.h"

ele_t get_ele_t(int ndim, int nnod/*, int neqn*/)		// only valid for equations implemented SO FAR
{														// later systems of eqns may violate these rules.
	ele_t type;
	if (ndim == 3)
	{
		type = iso3d;
	}
	else if (ndim == 2 && nnod == 4)
	{
		type = iso2d;
	}
	else if (ndim == 2 && nnod == 3)
	{
		type = tria2d;
	}
	else if (ndim == 1)
	{
		type = iso1d;
	}
	else
	{
		cout << "ERROR! Type not found!" << endl;
	}

	return type;
}

	// need to figure out how to resolve diff's between 2D and 3D functions when
	// they can't be determined by compiler -> templates? _3d appendage?

	// may want to come up with more generic algorithms
	// including find_all_faces_in_node_list
	// and include 2d, 3d relations.
	// also include reflective nodes (ie, local node #'s in ele across a face) etc.
int get_face( Tensor<bool, 1> & nodes )		// face using boolean matrix
{
	if ( nodes(0) && nodes(1) ) { return 0; }
	else if ( nodes(1) && nodes(2) ) { return 1; }
	else if ( nodes(2) && nodes(3) ) { return 2; }
	else if ( nodes(3) && nodes(0) ) { return 3; }
	else
	{
		cout << "error deterimining face using boolean node tensor: " << nodes << endl;
	//	cin.get();
		return -1;
	}
}

int get_face(int left, int right)
{
	if      ( left == 0 &&  right == 1 )	{ return 0; }	// Bottom face
	else if ( left == 1 &&  right == 2 )	{ return 1; }	// Right face
	else if ( left == 2 &&  right == 3 )	{ return 2; }	// Top face
	else if ((left == 0 && right == 3 ) 
		||  (right == 0 && left == 3) )		{ return 3; }	// Left face 
	else
	{
		cout << "error deterimining face using int node vector: " << left << " " << right << endl;
	//	cin.get();
		return -1;
	}
}

int get_face( std::vector<int> & nodes )	// face using node pair, any sequence.
{
	return get_face(nodes[0], nodes[1]);
}

int get_face( Tensor<int, 1> & nodes )	// face using node pair, any sequence.
{
	return get_face(nodes(0), nodes(1));
}

	// Returns the local node numbers of a face, given the face
Tensor<int, 1> get_nodes_from_face( int f)
{
	Tensor<int, 1> retval(2);
	retval(0) = f;
	retval(1) = f+1;
	if ( f+1 > 3 )
		retval(1) -= 4;

	return retval;
}

	// Returns the face number of the neighboring element at face f
int get_opposing_face( int f)
{
	f +=2;
	if ( f > 3 ) f -=4;				// shift face to the facing element -> might make into global function w/doc?
	return f;
}

/// 3D

	// Return nodes to form a face
	// may not be orderd by right-hand-rule for outward-pointing vector, verify later.
	// (I think it is, just not validated)
	// -> swapping ordering as shown to align with lists in MeshRefine. This orients the elements properly.
	//
Tensor<int, 1> get_nodes_from_face_3d( int f)
{
	Tensor<int, 1> face(4);

	if (f == 0)
	{
		face(0) = 3;		// try 1 0 3 2 
		face(1) = 2;
		face(2) = 1;
		face(3) = 0;
	}
	else if (f == 4)//(f == 1)
	{
		face(0) = 3;
		face(1) = 0;
		face(2) = 4;
		face(3) = 7;
	}
	else if (f == 1) //(f == 2)
	{
		face(0) = 0;
		face(1) = 1;
		face(2) = 5;
		face(3) = 4;	
	}
	else if (f == 2) // (f == 3)
	{
		face(0) = 1;
		face(1) = 2;
		face(2) = 6;
		face(3) = 5;
	}
	else if (f == 3) //(f == 4)
	{
		face(0) = 2;
		face(1) = 3;
		face(2) = 7;
		face(3) = 6;
	}
	else if (f == 5)
	{
/*		face(0) = 4;	// or 6 7 4 5
		face(1) = 5;	// need to "standardize" because we use 4 5 6 7 for bc's.
		face(2) = 6;
		face(3) = 7;
*/
		face(0) = 6;	// or 6 7 4 5
		face(1) = 7;
		face(2) = 4;
		face(3) = 5;

	}
	else
	{
		cout << "invalid face " << f << endl;
		cin.get();
	}

	return face;
}


int get_face_3d(int n0, int n1, int n2, int n3)
{
	if      ( n0 == 3 && n1 == 2 && n2 == 1 && n3 == 0 )	{ return 0; }	// Bottom face
	else if ( n0 == 3 && n1 == 0 && n2 == 4 && n3 == 7 )	{ return 4; }	// Front face
	else if ( n0 == 0 && n1 == 1 && n2 == 5 && n3 == 4 )	{ return 1; }	// Right face
	else if ( n0 == 1 && n1 == 2 && n2 == 6 && n3 == 5 )	{ return 2; }	// Back face
	else if ( n0 == 2 && n1 == 3 && n2 == 7 && n3 == 6 )	{ return 3; }	// Left face
	else if ( n0 == 4 && n1 == 5 && n2 == 6 && n3 == 7 )	{ return 5; }	// Top face
	else if ( n0 == 6 && n1 == 7 && n2 == 4 && n3 == 5 )	{ return 5; }	// Top face
	else
	{
		cout << "error deterimining face using int node vector: " << n0 << " " << n1 << " " << n2 << " " << n3 << endl;
		cin.get();
		return -1;
	}
}

int get_face_3d( std::vector<int> & nodes )	// face using node pair, any sequence.
{
//	std::sort(nodes.begin(), nodes.end(), sort_using_less_than);
	return get_face_3d(nodes[0], nodes[1], nodes[2], nodes[3]);
}

int get_face_3d( Tensor<int, 1> & nodes )	// face using node pair, any sequence.
{
//	vector<int> values;
//	values += nodes[0], nodes[1], nodes[2], nodes[3];
//	return get_face_3d(values);				// call vector version to get sort?
	return get_face_3d(nodes(0), nodes(1), nodes(2), nodes(3));
}


	// Returns the face number of the neighboring element at face f
int get_opposing_face_3d( int f)
{
	if ( f == 0)
		return 5;
	else if ( f == 1)
		return 3;
	else if ( f == 2)
		return 4;
	else if ( f == 3)
		return 1;
	else if ( f == 4)
		return 2;
	else if ( f == 5)
		return 0;
	else
	{
		cout << "No match in get_opposing_face_3d: " << f << endl;
		cin.get();
	}
}

	// get edge corresponding to nodelist
	// ni < nj
	// alternative would be declaring a nested array and iterating over it 
	// from 7 - 18 and use those values for ni, nj.
int get_edge_from_nodes(int ni, int nj)
{
	if ( ni > nj)
	{
		int larger = ni;
		ni = nj;
		nj = larger;
	}

	if (ni == 0 && nj == 1)
	{
		return 7;
	}
	else if (ni == 1 && nj == 2)
	{
		return 8;
	}
	else if (ni == 2 && nj == 3)
	{
		return 9;
	}	
	else if (ni == 0 && nj == 3)
	{
		return 10;
	}
	else if (ni == 0 && nj == 4)
	{
		return 11;
	}
	else if (ni == 1 && nj == 5)
	{
		return 12;
	}	
	else if (ni == 2 && nj == 6)
	{
		return 13;
	}	
	else if (ni == 3 && nj == 7)
	{
		return 14;
	}	
	else if (ni == 4 && nj == 5)
	{
		return 15;
	}	
	else if (ni == 5 && nj == 6)
	{
		return 16;
	}	
	else if (ni == 6 && nj == 7)
	{
		return 17;
	}	
	else if (ni == 4 && nj == 7)
	{
		return 18;
	}	
	else
	{
		cout << "no match in get_edge_from_nodes" << endl;
		cin.get();
		return 0;
	}
}

/*
int get_edge_from_nodes(int ni, int nj)
{
	if ( ni > nj)
	{
		int larger = ni;
		ni = nj;
		nj = larger;
	}

	if (ni == 2 && nj == 3)
	{
		return 7;
	}
	else if (ni == 3 && nj == 0)
	{
		return 8;
	}
	else if (ni == 0 && nj == 1)
	{
		return 9;
	}	
	else if (ni == 1 && nj == 2)
	{
		return 10;
	}
	else if (ni == 3 && nj == 7)
	{
		return 11;
	}
	else if (ni == 2 && nj == 6)
	{
		return 12;
	}	
	else if (ni == 1 && nj == 5)
	{
		return 13;
	}	
	else if (ni == 0 && nj == 4)
	{
		return 14;
	}	
	else if (ni == 6 && nj == 7)
	{
		return 15;
	}	
	else if (ni == 7 && nj == 4)
	{
		return 16;
	}	
	else if (ni == 4 && nj == 5)
	{
		return 17;
	}	
	else if (ni == 5 && nj == 6)
	{
		return 18;
	}	
	else
	{
		cout << "no match in get_edge_from_nodes" << endl;
		cin.get();
		return 0;
	}
}
*/


int get_face_tri( Tensor<int, 1> & nodes )		// face using boolean matrix
{
	if (	  nodes(0) == 0 && nodes(1) == 1 ) { return 0; }
	else if ( nodes(0) == 1 && nodes(1) == 2 ) { return 1; }
	else if ( nodes(0) == 2 && nodes(1) == 0 ) { return 2; }
	else
	{
		cout << "error deterimining face using node tensor: " << nodes << endl;
	//	cin.get();
		return -1;
	}
}

	// Returns the local node numbers of a face, given the face
Tensor<int, 1> get_nodes_from_face_tri( int f)
{
	Tensor<int, 1> retval(2);
	retval(0) = f;
	retval(1) = f+1;
	if ( f+1 > 2 )
		retval(1) -= 3;

	return retval;
}

// switcher

	// try this type of a system to clean up code ... 
	Tensor<int, 1> get_nodes_from_face( int face, ele_t ele_type )
	{
		if (ele_type == iso2d)
		{
			return get_nodes_from_face( face);
		}
		else if ( ele_type == iso3d)
		{
			return get_nodes_from_face_3d( face);
		
		}
		else if (ele_type == tria2d)
		{
			return get_nodes_from_face_tri( face);
		}
		else if (ele_type == iso1d)
		{
			Tensor<int,1> retval(1);
			retval(0) = face;
			return retval;
		}
	}

	int get_face_1d( Tensor<int, 1> &n)
	{
		return n(0);
	}

	int get_face( Tensor<int, 1> & n, ele_t ele_type )
	{
		cout << "type: " << ele_type << " size: " << n.imax() << endl;
		if (ele_type == iso2d)
		{
			return get_face( n);
		}
		else if ( ele_type == iso3d)
		{
			return get_face_3d( n);
		}
		else if (ele_type == tria2d)
		{
			return get_face_tri( n);
		}
		else if (ele_type == iso1d)
		{
			return get_face_1d( n);
		}
	}

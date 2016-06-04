#pragma once

#include "FDV/Node.h"
//#include "FDV/BC/BoundaryCondition.h"
// holds face information for boundary conditions
struct Face
{
	~Face()
	{
//		std::cout << "~Face(): bc " << bc << " n.size() " << n.imax << std::endl;
	}
	Face() : n(0) {};
	Face (int nbnod)
		: n(nbnod)
	{
		bc = 0;
	}
	int face;								// face the bc applies to (redundant with node list i guess)
	Tensor<int, 1> n;
	int bc;									// face boundary condition
};


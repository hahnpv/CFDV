#pragma once

/// base class for GMRES solvers
/// -> MPI
/// -> unparallelized
/// -> CUDA

struct GMRESbase
{
	/// Constructor, intializes dimensions of arrays
	GMRESbase()
	{}

	/// Solves the system of equations, given the element correlation and node data
	/// @param elements element data structure, provides input data
	virtual void iterate( std::vector< Element *> & elements )
	{}

	/// Update the nodes with the residuals
	/// @param nodes NodeVar data structure, for updating results
	virtual void update( std::vector< Node *> & nodes )
	{}
};


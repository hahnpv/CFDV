#pragma once

#include "Point.h"


struct Line : public Tensor<double, 2>
{
	Line( Point & left, Point & right)
		: Tensor<double, 2>(2,2)					// 2, ndim
	{
		(*this)(0, 0) = left(0);
		(*this)(0, 1) = left(1);

		(*this)(1, 0) = right(0);
		(*this)(1, 1) = right(1);
	}

	Line( Point & left, Line & right)
		: Tensor<double, 2>(right.imax()+1, 2)
	{
		(*this)(0, 0) = left(0);
		(*this)(0, 1) = left(1);

		for (int i = 1; i < right.imax() + 1; i++)
		{
			(*this)(i, 0) = right( i - 1, 0);
			(*this)(i, 1) = right( i - 1, 1);
		}
	}

	Line( Line & left, Point & right)
		: Tensor<double, 2>(right.imax()+1, 2)
	{
		for (int i = 0; i < right.imax(); i++)
		{
			(*this)(i, 0) = left( i, 0);
			(*this)(i, 1) = left( i, 1);
		}

		(*this)(right.imax(), 0) = right(0);
		(*this)(right.imax(), 1) = right(1);
	}

	// need to think about operation w/VVF, do we subsume the functionality here 
	// or do we re-construct VVF to work with Line?

private:
};

		// Merges two points into a line
Line operator+(Point left, Point right)	
{
	return Line(left, right);
}

		// Inserts a point into a line
Line operator+(Point left, Line right)
{
	return Line(left, right);
}

		// Appends a point to a line
Line operator+(Line left, Point right)
{
	return Line(left, right);
}
#pragma once

#include "../Utility/Tensor.h"

	// Convenience function for a Tensor which represents a 2D or 3D point in space
struct Point : public Tensor<double, 1>
{
	Point(double x, double y)
		: Tensor<double, 1>(2)
	{
		(*this)(0) = x;
		(*this)(1) = y;
	}

	Point(double x, double y, double z)
		: Tensor<double, 1>(3)
	{
		(*this)(0) = x;
		(*this)(1) = y;
		(*this)(2) = z;
	}
};
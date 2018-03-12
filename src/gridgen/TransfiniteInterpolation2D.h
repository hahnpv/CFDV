#pragma once

	// Could submit xi[], eta[] interpolation functions from [0,1]
template< class T>
Tensor<T, 3> TransfiniteSurface(Tensor<T, 2> & bottom, Tensor<T, 2> & top, Tensor<T, 2> & left, Tensor<T, 2> & right, int delta_i, int delta_j)
{
	VectorValuedFunction Fb(bottom, 0);			// looks ok at a glance
	VectorValuedFunction Ft(top,    0);			// looks ok at a glance
	VectorValuedFunction Fl(left,   1);			// looks ok at a glance
	VectorValuedFunction Fr(right,  1);			// looks ok at a glance

		// **must ensure** connectivity between lines at these points.
	Point F00(bottom(0,0), bottom(0,1));		// bottom left
	Point F01(right(0,0), right(0,1));			// bottom right
	Point F10(top(0,0), top(0,1));				// top left
	Point F11(right(right.imax()-1,0), right(right.imax()-1,1));		// top right

	Tensor<double, 3> U(delta_i+1,delta_j+1,2);		// domain

	for (int i = 0; i <= delta_i; i++)							// i == xi
	{
		double xi = (double)i / (double)delta_i;

		for (int j = 0; j <= delta_j; j++)						// j == eta
		{
			double eta = (double)j / (double)delta_j;

			for (int k = 0; k < 2; k++)							// (x,y)
			{
				U(i, j, k) = (1. - xi) * Fl(eta, k) + xi * Fr(eta, k) + (1 - eta) * Fb(xi, k) + eta * Ft(xi, k)
					- ( (1.-xi)*(1.-eta)*F00(k) + (1.-xi)*eta*F10(k) + xi * (1.-eta)*F01(k) + xi * eta * F11(k) );

			}
		}
	}

	return U;
}

#pragma once

	/// Clears the transient data in each element
template<class T> struct AxisymmetricFlow : public unary_function<T, void>
{
	AxisymmetricFlow(int neqn, int nnod)
		: neqn(neqn), nnod(nnod)
	{
	}
	void operator() (T e) 
	{
		// According to Chung (verify), a 2D axisymmetric flow about x=0 is calculated by
		// 2 * pi * r * integral(xi, eta),
		// so should in theory be able to postmultiply a 2D FDV galerkin calculation by 2*PI (with r in or out of integral?).

		double pi = 3.1415926535857;

		double r = 0;
		for (int i = 0; i < 4; i++)
		{
			r += (1.0/4.0) * e->node[i]->x(1);
		}

		for (int i = 0; i < neqn * nnod; i++)
		{
			for (int j = 0; j < neqn * nnod; j++)
			{
				e->R(i, j) *= 2 * pi;// * r;
			}

			e->rhs(i) *= 2 * pi;// * r;
		}

	}
	int neqn;
	int nnod;
};


#pragma once

	/// Extracts primitive variables from the conservation variable
template<class T> struct NodeUnpack //: public unary_function<T, void>
{
	NodeUnpack()
	{
	}
	void operator() (T node)
	{
		Thermo & thermo = Thermo::Instance();

		clear(node);

		node->rho  = node->U(0);

		double vsqr = 0;
		for (unsigned int i=0; i<node->x.imax(); i++)
		{
			node->v(i) = node->U(1+i) / node->U(0);
			vsqr += (node->v(i) * node->v(i));				// could use dot() out of loop
		}

		node->E = node->U(node->x.imax()+1) / node->U(0);
		node->e = node->E - 0.5 * vsqr;

		node->T = node->e / thermo.Cv;
		node->p = node->rho * (thermo.gamma-1.0) * node->e;
	}

		// zero out all flow transients
	void clear(T node)
	{
		node->rho = 0;

		node->v.clear();

		node->p = 0;
		node->T = 0;

		node->E = 0;
		node->e = 0;
	}

};


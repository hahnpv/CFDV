#pragma once

struct NavierStokes
{
	// 2D Convective Flux
	virtual void calc_F(Tensor<double, 1> & U, double p, Tensor<double, 2> & matrix) {};

	// 2D Viscous Flux
	virtual void calc_G(Tau & tau, Tensor<double, 1> & v, Tensor<double, 1> & dT, double & HK, double & Re, Tensor<double, 2> & matrix) {};

	/// 2D Source Term, not incorportated
	virtual void calc_B(Node * n, Tensor<double,1> & matrix) {};

	// 2D Flux Jacobians

	// dF/dU
	virtual void calc_a(Flow & f, Tensor<double,3> & matrix) {};

	// dG/dU
	virtual void calc_b(Flow & flow, Tensor<double, 3> & matrix) {};

	// dG/dU,j		
	virtual void calc_c(Flow & flow, Tensor<double, 4> & matrix) {};

	// dB/dU
	virtual void calc_d(Tensor<double, 2> & matrix) {};
};


#pragma once

#include <iostream>
#include <fstream>
#include "RMSBase.h"
using namespace std;

	/// monitor RMS Error of CFD solution, MPI parallelized
template<class T> struct MPI_RMSError : public RMSBase<T> // public unary_function<T, void>
{
	MPI_RMSError(std::string filename, double t, int iter, int neqn)
		: filename( filename), t(t), iter(iter), dU(neqn), dU0(neqn), neqn(neqn), Urms(neqn)
	{
		dU.clear();
		dU0.clear();
	}

	typedef void result_type;					// result_type for tr1 wrapper.
	result_type operator() (T n) 
	{
		for (int j=0; j < neqn; j++)
		{
			dU(j)  += n->U(j)  *  n->U(j);
			dU0(j) += n->U0(j) * n->U0(j);
		}
	}
	~MPI_RMSError()
	{
		Tensor<double, 1> Urms(neqn);
		Urms.clear();

		Tensor<double, 1> dUlocal(neqn);
		Tensor<double, 1> dU0local(neqn);
		for (int i = 0; i < neqn; i++)
		{
			dUlocal(i)  = dU(i);
			dU0local(i) = dU0(i);
		}

		for (int i = 0; i < neqn; i++)
		{
			MPI::COMM_WORLD.Allreduce(&dUlocal(i), &dU(i), 1, MPI_DOUBLE, MPI_SUM);
			MPI::COMM_WORLD.Allreduce(&dU0local(i), &dU0(i), 1, MPI_DOUBLE, MPI_SUM);
		}

		for (int i = 0; i < neqn; i++)
		{
			Urms(i) = sqrt( dU(i) / dU0(i) ) - 1.0;
		}

		if ( MPI::COMM_WORLD.Get_rank() == 0)
		{
			ofstream fout(filename.c_str(), ios::app);

			fout << iter << "\t";
			cout << "Residual: ";

			for (int i=0; i < neqn; i++)
			{
				fout << Urms(i) << " ";
				cout << Urms(i) << " ";
			}

			fout << endl;
			cout << endl;

			fout.close();
		}
	}

	Tensor<double, 1> Urms;
	Tensor<double, 1> dU;
	Tensor<double, 1> dU0;

	double t;
	int iter;
	int neqn;

	std::string filename;
};


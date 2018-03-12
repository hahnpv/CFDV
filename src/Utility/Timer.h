#pragma once

// a timer for CFD iterations

struct IterationTimer
{
	IterationTimer(int rank)
		: running_time(0),
		  rank(rank)
	{
		if (rank == 0)
		{
			ofstream fout("time.txt",ios::out);
			fout << "Iteration, Run_Time, Time_Elapsed, Sim_Time, Step_dt" << endl;
			fout.close();
		}
	}

	void reset()
	{
		timr.restart();
	}

	void check(int iter, double t, double dt)
	{
		if (rank == 0)
		{
			running_time += timr.elapsed();
			ofstream fout("time.txt",ios::app);
			fout << iter << " " << running_time << " " << timr.elapsed() << " " << t << " " << dt << " " << endl;
			fout.close();
		}
	}

	int rank;
	double running_time;
	boost::timer timr;
};
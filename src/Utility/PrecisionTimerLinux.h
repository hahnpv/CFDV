#include <sys/time.h>
#include <iostream>
#include <vector>
using namespace std;

struct LinuxPerformanceTimer
{
	void start()
	{
		time0 = get_time();
	}

	void stop()
	{
		time1 = get_time();
	}

	double read()
	{
		return (time1 - time0) * 1E-9;
	}
	
private:
	unsigned long long get_time()
	{	
		#define rdtsc(low,high) \
			__asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high))

		unsigned long low, high;
		rdtsc(low, high);
		return ((unsigned long long)high << 32) | low;
	}

	unsigned long long time0;
	unsigned long long time1;
};


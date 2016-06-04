#pragma once

// Windows-specific precision timer
// need to add in Linux code and #ifdef block it by machine type



#ifdef __GNUC__
	#include <sys/time.h>
#else
	#include <tchar.h>
	#include <windows.h>
#endif

/*
	Use:

//		PrecisionTimer timer;
//		timer.start();
//		do something		
//		timer.stop();
//		cout << "execution time: " << timer.read() << endl;

*/

#ifdef __GNUC__
	struct PT
	{
		LinuxPT() {};
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
#else
	struct PT
	{
	public:
		PT()
		{
			QueryPerformanceFrequency(&ticks_per_second);
		}

		void start()
		{
			QueryPerformanceCounter(&tick1);	
		}

		void stop()
		{
			QueryPerformanceCounter(&tick2);	
		}

			// Returns the time in seconds
		double read()
		{
			return (double)(tick2.QuadPart - tick1.QuadPart) / (double)ticks_per_second.QuadPart;
			return 0;
		}

	private:
		LARGE_INTEGER ticks_per_second;
		LARGE_INTEGER tick1;
		LARGE_INTEGER tick2;
	};
#endif

class PrecisionTimer
{
public:
	PrecisionTimer() {};

	void start()
	{
//		QueryPerformanceCounter(&tick1);	
		timer.start();
	}

	void stop()
	{
//		QueryPerformanceCounter(&tick2);	
		timer.stop();
	}

		// Returns the time in seconds
	double read()
	{
//		return (double)(tick2.QuadPart - tick1.QuadPart) / (double)ticks_per_second.QuadPart;
//		return 0;
		return timer.read();
	}

private:
	PT timer;
};


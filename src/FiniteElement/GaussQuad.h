// gaussian quadrature, 1D

#pragma once

//#include "TestFunctions1D.h"

// could use Boost to assign w, zetas in one line maybe?

template <class TClass> class GaussQuad
{
public:
	GaussQuad(int N)
		: N(N),
		  w(N),
		  zeta(N)
	{ 
		get_gauss();
	};

	~GaussQuad()
	{}

	/* // Loop Unrolling Example
	template< int i >
	class LOOP{
	  public:
		static inline void EXEC(){
		  cout << "A-" << i << " ";
				LOOP< i-1 >::EXEC();
		   cout << "B-" << i << " ";
		}
	};

	class LOOP< 0 >{
	  public:
		static inline void EXEC(){
		  cout << "A-" << i;
		  cout << "\n"; 
		   cout << "B-" << i;
		}
	};

	template<> // need to know N at compile time!
	class
	{
		public:
		static inline void EXEC(w, zeta)
		{
			callback ptr
		};
	};
	*/

		/// One Dimensional Quadrature
	virtual void one(TClass* pt2Object, void(TClass::*fpt)(Tensor<double,1> &,Tensor<double,1> &))
	{
		Tensor<double,1> _w1( 1);
		Tensor<double,1> _zeta1( 1);
		for (unsigned int i=0; i < N; i++)
		{
			_w1(0) = w(i);
			_zeta1(0) = zeta(i);

			(*pt2Object.*fpt)(_w1, _zeta1);			// execute member function
		}
	};  

		/// Two Dimensional Quadrature
	virtual void two(TClass* pt2Object, void(TClass::*fpt)(Tensor<double,1> &,Tensor<double,1> &))
	{
		Tensor<double,1> _w2( 2);
		Tensor<double,1> _zeta2( 2);
		for (unsigned int i=0; i < N; i++)
		{
			for (unsigned int j=0; j < N; j++)
			{
				_w2(0) = w(i);
				_w2(1) = w(j);

				_zeta2(0) = zeta(i);
				_zeta2(1) = zeta(j);

				(*pt2Object.*fpt)(_w2, _zeta2);			// execute member function
			}
		}
	};

		/// Three Dimensional Quadrature
	virtual void three(TClass* pt2Object, void(TClass::*fpt)(Tensor<double,1> &,Tensor<double,1> &))
	{
		Tensor<double,1> _w3( 3);
		Tensor<double,1> _zeta3( 3);

		for (unsigned int i=0; i < N; i++)
		{
			for (unsigned int j=0; j < N; j++)
			{
				for (unsigned int k=0; k < N; k++)
				{				
					_w3(0) = w(i);
					_w3(1) = w(j);
					_w3(2) = w(k);

					_zeta3(0) = zeta(i);
					_zeta3(1) = zeta(j);
					_zeta3(2) = zeta(k);

					(*pt2Object.*fpt)(_w3, _zeta3);			// execute member function
				}
			}
		}
	};  


		/// Generate w, zeta for a given N
	void get_gauss()
	{
		double x1 = -1.0;		// area of integration
		double x2 =  1.0;	
		double eps = 3.E-14;

		double m=(N+1)/2;				//  1.5
		double xm=0.5*(x2+x1);		//	0
		double xl=0.5*(x2-x1);		//  1

		Tensor<double,1> p(4);
		for (int i=1; i <= m; i++)
		{
			double z = cos(3.141592654*(i-0.25)/(N+0.5));

			double z1 = z + 10;		// initial error to avoid satisfying while loop
			p(0) = 0;
			while(abs(z-z1) > eps)
			{
				p(1)=1.0;
				p(2)=0.0;
				for (unsigned int j=1; j <= N; j++)
				{
					p(3)=p(2);
					p(2)=p(1);
					p(1)=((2.0*j-1.0)*z*p(2)-(j-1.)*p(3))/j;
				}
				p(0)=N*(z*p(1)-p(2))/(z*z-1.0);
				z1=z;
				z=z1-p(1)/p(0);
		   }
		   zeta(i-1)=xm-xl*z;
		   zeta((unsigned int)N+1-i-1)=xm+xl*z;						// -1 c++ count
		   w(i-1)=2.0*xl/((1.0-z*z)*p(0)*p(0));
		   w((unsigned int)N+1-i-1)=w(i-1);							// -1 c++ count
		}
		for(unsigned int i=0; i < N; i++)
		{
			if(abs(zeta(i)) < eps) 
				zeta(i)=0.0;
			if(abs(w(i)) < eps) 
				w(i)=0.0;
		}
	}

private:
	unsigned int N;
	Tensor<double, 1> w;
	Tensor<double, 1> zeta;
};


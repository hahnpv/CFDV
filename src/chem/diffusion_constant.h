#pragma once


// Diffusion constant for rigid spheres of unequal mass (m_a, m_b) and diameter (d_a, d_b)

// T = molecule

template < class molecule>
double Dab( molecule a, molecule b, double T, double p)
{
	double d_a = a.diameter;
	double d_b = b.diameter;

	double m_a = a.mass;
	double m_b = b.mass;

	double pi = 3.14159265;
	double k  = 1.3806503 * 10E23;			// Boltzmann constant, do we need to nondimensionalize?

	return (2./3.) * sqrt(pow(k, 3) / pow(pi,3)) * sqrt( ( 1 / (2. * m_a)) + ( 1 / (2. * m_b)) ) * ( sqrt(T) / ( p * pow( (d_a + d_b) / 2.0, 2 ) );
}
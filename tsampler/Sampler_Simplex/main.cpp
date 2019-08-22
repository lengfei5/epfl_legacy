#include "Sampler.h"
#include "Test.h"

// g++ -m32 main.cpp Test.cpp -o Sampler -lgsl -I/sw/include

int main (int argc, char * const argv[])
{
	const size_t N_it  = 10000; // number of iteration
	const size_t N_var = 3;     // number of variable
	double * alpha = new double[N_var]; // alpha values for prior (dirichlet)
	double * P_ini = new double[N_var]; // initial values of the parameters
	
	alpha[0] = 5.;
	alpha[1] = 5.;
	alpha[2] = 5.;
	
	P_ini[0] = 0.2;
	P_ini[1] = 0.3;
	P_ini[2] = 0.5;
	
	// test object, owns the method "double compute_likelihood(const double* const P)"
	// where P is a vector of parameters of size N_var
	test t1;
	
	// construction of the sampler, takes as an argument an pointer on an object 
	// which owns a method "double compute_likelihood(const double* const P)"
	sampler<test> mysampler(&t1,3,5.,alpha);
	// Sampling process
	mysampler.run_sampling(P_ini, N_it);
	// Print values in a file, i.e. parameters, liklihood and if the move was accepted
	// or not at each iteration
	mysampler.print_result("test.dat");
	
	delete [] alpha;
	delete [] P_ini;
}

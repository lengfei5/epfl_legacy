/*
 *  Sampler.h
 *  Sampler_Simplex
 *
 *  Created by Benjamin Zoller on 01.11.12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 *  Class template sampler
 *
 */

#ifndef _sampler_h_included_
#define _sampler_h_included_

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../gene.hh"
#include "../output.hh"

using namespace std;

template <typename T> 
class sampler
{
private:
	size_t N_var;
	size_t N_it;
	double scale;
	bool paired;
	
	double *  alpha_prior;
	double ** Values;
	
	T * myobject;
	
	gsl_rng * r;
	
	double Proposal_Ratio(const double * const P_prev, const double * const P_next);
	double Prior_Ratio(const double * const P_prev, const double * const P_next);
	void Move_Proposal(const double * const P_prev, double * const P_next);
	size_t Move_Acceptance(const double * const P_prev, const double * const P_next, const double Likelihood_prev, const double Likelihood_next);
	
public:
	sampler()
	{
		N_var = 0;
		N_it  = 0;
		scale = 0.;
		
		alpha_prior = NULL;
		Values      = NULL;
		
		myobject = NULL;
		r = NULL;
	}
	
	sampler(T* myobject_, const size_t N_var_, const double scale_, const double * const alpha_, bool _p)
	{
		paired = _p;
		N_var = N_var_;
		N_it  = 0;
		scale = scale_;
		
		alpha_prior = new double[N_var];
		for (size_t k = 0; k < N_var; ++k) alpha_prior[k] = alpha_[k];
		
		// store evolution of parameters, likelihood and acceptance ratio //
		
		Values = new double*[N_var+2];
		
		// store pointer on object for Likelihood computation //
		
		myobject = myobject_;
		
		// seed for random generator //
		
		int seed = time(0);
		r = gsl_rng_alloc (gsl_rng_mt19937);
		gsl_rng_set (r, seed);
	}
	
	~sampler()
	{
		
		delete [] alpha_prior;
		
		if (N_it > 0)
		{
			for (size_t k = 0; k < N_var+2; ++k) delete [] Values[k];
		}
		delete [] Values;
		gsl_rng_free (r);
	}
	
	void run_sampling(const double * const P_ini, const size_t N_it);
	void print_result(string fname);
	void print_summary(int burnin, gene &g);
	void get_summary(int burnin, output_one_condition &out);
};

template <typename T>
double sampler<T>::Proposal_Ratio(const double * const P_prev, const double * const P_next)
{
	double * alpha_next = new double[N_var];
	for (size_t i = 0; i < N_var; ++i) alpha_next[i] = 10e-3 + scale*P_next[i];
	
	double * alpha_prev = new double[N_var];
	for (size_t i = 0; i < N_var; ++i) alpha_prev[i] = 10e-3 + scale*P_prev[i];

	const double Proposal_prev = gsl_ran_dirichlet_pdf (N_var, alpha_next, P_prev);
	const double Proposal_next = gsl_ran_dirichlet_pdf (N_var, alpha_prev, P_next);
	
	delete [] alpha_next;
	delete [] alpha_prev;
	
	return Proposal_prev/Proposal_next;
}

template <typename T>
double sampler<T>::Prior_Ratio(const double * const P_prev, const double * const P_next)
{
	const double Prior_prev = gsl_ran_dirichlet_pdf (N_var, alpha_prior, P_prev);
	const double Prior_next = gsl_ran_dirichlet_pdf (N_var, alpha_prior, P_next);
	
	return Prior_next/Prior_prev;
}

template <typename T>
void sampler<T>::Move_Proposal(const double * const P_prev, double * const P_next)
{
	double * alpha_next = new double[N_var];
	for (size_t i = 0; i < N_var; ++i) alpha_next[i] = 10e-3 + scale*P_prev[i];
	
	gsl_ran_dirichlet (r, N_var, alpha_next, P_next);
	
	delete [] alpha_next;
}

template <typename T>
size_t sampler<T>::Move_Acceptance(const double * const P_prev, const double * const P_next, const double Likelihood_prev, const double Likelihood_next)
{
	const double Likelihood_Ratio = exp(Likelihood_next - Likelihood_prev); // log likelihood
	const double Acceptance_Prob = Likelihood_Ratio * Prior_Ratio(P_prev, P_next) * Proposal_Ratio(P_prev, P_next);
	
	size_t accepted = 0;
	if (gsl_rng_uniform (r) < Acceptance_Prob) accepted = 1;
	
	return accepted;	
}

template <typename T>
void sampler<T>::run_sampling(const double * const P_ini, const size_t N_it_)
{
	// Initialization //
	
	N_it = N_it_;
	
	for (size_t i = 0; i < N_var+2; ++i)
	{
		Values[i] = new double[N_it];
	}
	
	double * P_prev = new double[N_var];
	double * P_next = new double[N_var];
	
	for (size_t i = 0; i < N_var; ++i)
	{
		P_prev[i]    = P_ini[i];
		Values[i][0] = P_ini[i];
	}
	
	double Likelihood_prev;
	double Likelihood_next = 0.;
	if (paired)
	{ 
	Likelihood_prev = myobject->compute_likelihood_paired(P_prev);
	}
	
	else
	{
	Likelihood_prev = myobject->compute_likelihood_single(P_prev);
	}

	Values[N_var][0]   = Likelihood_prev;
	Values[N_var+1][0] = 1.;
	
	// Begin sampling //
	
	size_t accepted = 0;
	
	for (size_t it = 1; it < N_it; ++it)
	{
		Move_Proposal(P_prev, P_next);
		
		//cerr << it << " " << P_next[0] << " " << P_next[1] << " " << P_next[2] << endl;
		if (paired)
		{
		Likelihood_next = myobject->compute_likelihood_paired(P_next);
		}
		else
		{
		 Likelihood_next = myobject->compute_likelihood_single(P_next);
		}
		accepted = Move_Acceptance(P_prev, P_next, Likelihood_prev, Likelihood_next);
		if (accepted)
		{
			for (size_t k = 0; k < N_var; ++k)
			{
				P_prev[k]     = P_next[k];
				Values[k][it] = P_next[k];
			}
			
			Values[N_var][it]   = Likelihood_next;
			Values[N_var+1][it] = accepted;
			Likelihood_prev = Likelihood_next;
		}
		else
		{
			for (size_t k = 0; k < N_var; ++k)
			{
				Values[k][it] = P_prev[k];
			}
			
			Values[N_var][it]   = Likelihood_prev;
			Values[N_var+1][it] = accepted;
		}
	}
	
	delete [] P_prev;
	delete [] P_next;
}

template <typename T>
void sampler<T>::print_result(string fname)
{
	ofstream ofs(fname.c_str());
	
	for (size_t it = 0; it < N_it; ++it)
	{
		for (size_t k = 0; k < N_var+1; ++k)
		{
			char out[128];
			// vector(out, "%.3e\t", Values[k][it]);
			// ofs << Values[k][it] << " ";
			ofs << out;
		}
		char out[128];
		sprintf(out, "%.f", Values[N_var+1][it]);
		ofs << out;
		ofs << endl;
	}
	
	ofs.close();
}

template <typename T>
void sampler<T>::print_summary(int burnin, gene &g)
{
	vector<double> sum(N_var+2,0.);
	vector<double> sum2(N_var+2,0.);
	unsigned long int n=0;
	
	for (size_t it = burnin; it < N_it; ++it)
	{
		for (size_t k = 0; k < N_var+2; ++k)
		{
			sum[k]  += Values[k][it];
			sum2[k] += Values[k][it]*Values[k][it];
		}
		n++;
	}

	// printf("%s\t", prefix);
	for (size_t k = 0; k < N_var; ++k)
	{
		printf("%s\t%s\t%d\t%.6f\t%.6f\t%.3e\t%.2f\n", g.name.c_str(), g[k].name.c_str(), g.nreads(), sum[k]/n, sqrt(sum2[k]/n - (sum[k]/n)*(sum[k]/n)),  sum[N_var]/n, sum[N_var+1]/n);
	}
	// int k=N_var;
	// printf("%.3e\t", sum[k]/n);
	// k++;
	// printf("%.3f\n", sum[k]/n);
}

template <typename T>
void sampler<T>::get_summary(int burnin, output_one_condition &out)
{
	vector<double> sum(N_var+2,0.);
	vector<double> sum2(N_var+2,0.);
	unsigned long int n=0;
	
	for (size_t it = burnin; it < N_it; ++it)
	{
		for (size_t k = 0; k < N_var+2; ++k)
		{
			sum[k]  += Values[k][it];
			sum2[k] += Values[k][it]*Values[k][it];
		}
		n++;
	}

	for (size_t k = 0; k < N_var; ++k)
	{
		out[k].mean = sum[k]/n;
		out[k].sd   = sqrt(sum2[k]/n - (sum[k]/n)*(sum[k]/n));
	}
	out.acceptance_ratio    = sum[N_var+1]/n;
	out.mean_log_likelihood = sum[N_var]/n;
}

#endif
//test

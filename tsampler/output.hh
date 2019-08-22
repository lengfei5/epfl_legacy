#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>

#include "gene.hh"
#include "output.hh"
#include "bamfiles.hh"
using namespace std;

typedef struct {
	double mean;
	double sd;
} summary;

class output_one_condition : public vector<summary>{
	//the different trancripts
public:

	output_one_condition(int _nt, int _nreads):vector<summary>(_nt), nreads(_nreads){};
	
	int nreads;
	//string condition_name;
	double acceptance_ratio;
	double mean_log_likelihood;
};

class output_all : public vector <output_one_condition>
	//all conditions
{
public:
	output_all(gene &_g):g(_g){};
	string genename;
	void print(unsigned int ng){
		
		//print header if first time;
		if(ng==0)
		{
		//	cout << "gene" << "\t" << "transcript" << "\t";
			for(unsigned int n=0; n< size(); ++n)//n loops over conditions
			{
				//printf("reads:%s\tmean:%s\tsd:%s\tar:%s\t", at(n).condition_name.c_str(),at(n).condition_name.c_str(),at(n).condition_name.c_str(),at(n).condition_name.c_str());
			}
			printf("\n");
			// printf("%.3e\t%.3f\n", at(size()-1).mean_log_likelihood, at(size()-1).acceptance_ratio);
		}
	
		for(unsigned int m=0; m<at(0).size(); ++m)//m loops over transcripts
		{	
			cout << g.name << "\t" << g[m].name << "\t";
			for(unsigned int n=0; n<size(); ++n)//n loops over conditions
			{
				printf("%d\t%.3f\t%.3f\t%.3f\t", at(n).nreads, at(n).at(m).mean, at(n).at(m).sd, at(n).acceptance_ratio);
			}
			printf("\n");
			// printf("%.3e\t%.3f\n", at(size()-1).mean_log_likelihood, at(size()-1).acceptance_ratio);
		}
	};
	
	void print_only_total(unsigned int ng)
	{	
		//print header if first time;
		if(ng==0)
		{
			cout << "gene" << "\t";
			for(unsigned int n = 0; n<size(); ++n)//n loops over conditions
			{
				// printf("reads:%s\t", bams[n].name);
			}
			printf("\n");
			// printf("%.3e\t%.3f\n", at(size()-1).mean_log_likelihood, at(size()-1).acceptance_ratio);
		}
	
		cout << g.name << "\t";
		for(unsigned int n  = 0; n<size(); ++n)//n loops over conditions
		{
			printf("%d\t", at(n).nreads);
		}
		printf("\n");
	};


private:
	gene &g; //CHECK THIS
	// bamfiles &bams;
};

#endif

#ifndef GENES_HH
#define GENES_HH

#include <vector>
#include <map>
#include <string>
#include "bamfiles.hh"
using namespace std;

typedef unsigned long long int transcript_index;

class transcript
{
public:
	transcript(string _name, int _i, int _l, double _w):name(_name), index(_i), len(_l),weight(_w){};
	string name;
	unsigned int index;
	unsigned int len;
	double prob;
	double weight;
};

class hit
{
public:

	hit(transcript* _t, unsigned int _p,unsigned int _is):trans(_t),pos(_p), isize(_is){};

	transcript* trans;
	unsigned int pos; // zero-based position on transcript
	unsigned int isize;
};

unsigned int number_of_keys(multimap<string, hit> &mmap);


class gene: public vector<transcript>
{
public:
	gene(string _name, vector<string> &tnames, int _rl,vector<double> frag);
	void import_reads(vector<bamfile> &bam);
	double compute_likelihood_paired(double *);
	double compute_likelihood_single(double *);
	void print();
	void setprobs(double *p);
	
	void classify_single_reads();
	void classify_paired_reads();
	unsigned int nreads(){return number_of_keys(mapped_reads);}

	string name;
	int readlen;
	
private:
	multimap<string, hit> mapped_reads;
	map<transcript_index,unsigned int> number_per_class;
	vector<double> frag_distribution;
	vector< vector<hit> > reads;
};

#endif

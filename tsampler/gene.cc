#include "gene.hh"

#include <assert.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "samtools/sam.h"
#include "tokens.hh"
#include "bamfiles.hh"
typedef struct {
    transcript *t;
    multimap<string, hit> *m;
} bamdata;

unsigned int number_of_keys(multimap<string, hit> &mmap)
{
	multimap<string, hit>::iterator m_it, s_it;
	unsigned int num=0;
	for(m_it = mmap.begin(); m_it!=mmap.end(); m_it=s_it)
	{
	    string theKey = (*m_it).first;
		// cout << theKey << endl;
	    pair<multimap<string, hit>::iterator,multimap<string, hit>::iterator> keyRange = mmap.equal_range(theKey);

		s_it = keyRange.second;
		num++;
	}
	return num;
}

gene::gene(string _name, vector<string> &tnames, int _rl, vector<double> frag):name(_name), readlen(_rl)
{	
	frag_distribution = frag;
	unsigned int index=0;
	for(vector<string>::iterator it=tnames.begin(); it !=tnames.end(); ++it)
	{
		vector<string>tmp;
		TokenizeV(*it, tmp, "|");
		unsigned int length = atoi(tmp[1].c_str());
		
		double w = 0.;
			
		if(frag_distribution.size() != 0 && length > 101){
		
		for(unsigned int i = 1; i <= (length - readlen + 1) ;  ++i){
			for(unsigned int j = readlen ; j <= (length - i + 1); ++j){
				
				if( j < 800){
					w += frag_distribution[j - readlen];
			
				}
			}
		}
}
	
                 w = 0;
		push_back(transcript(tmp[0], index, length, w));
		index++;
	}
};

void gene::setprobs(double *p)
{
	for(unsigned int n=0; n<size(); ++n)
	{
		at(n).prob = p[n];
	}
};


void gene::classify_paired_reads()
{
	reads.clear();
	int index = 0;	
	multimap<string, hit>::iterator m_it, s_it;

	for(m_it = mapped_reads.begin(); m_it!=mapped_reads.end(); m_it=s_it)
	{
		string theKey = (*m_it).first;
		// cout << theKey << endl;
		pair<multimap<string, hit>::iterator,multimap<string, hit>::iterator> keyRange = mapped_reads.equal_range(theKey);
		vector<hit> temp_hit;

		for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
		{
			index++;
			temp_hit.push_back(s_it->second);

		}
		reads.push_back(temp_hit);
	}

//cout << index << endl;
};

void gene::classify_single_reads()
{
		
	// cout << sizeof(transcript_index) << endl;
	cout << size() << endl;
	cout << sizeof(transcript_index) << endl;
	if(size()>8*sizeof(transcript_index)) {cout << "number of transcripts larger than available bits in index\n"; exit(1);}
		
	number_per_class.clear();
	
	multimap<string, hit>::iterator m_it, s_it;
	for(m_it = mapped_reads.begin(); m_it!=mapped_reads.end(); m_it=s_it)
	{
	    string theKey = (*m_it).first;
		// cout << theKey << endl;
	    pair<multimap<string, hit>::iterator,multimap<string, hit>::iterator> keyRange = mapped_reads.equal_range(theKey);

		int class_index=0;
		for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
		{
			transcript *trans = s_it->second.trans;
			class_index |= (1 << trans->index);
			// cout << "CI:" << class_index << endl;
		}
		number_per_class[class_index]++;
	}
	//check sum
	unsigned int total=0;
	for(map<transcript_index, unsigned int>::iterator it=number_per_class.begin(); it!=number_per_class.end(); ++it)
	{
		// cout << it->first <<":" << it->second << endl;
		total += it -> second;
	}

//	cout << total << endl;

	// cout << total << endl;

	// cout << mapped_reads.size() << endl;
	if(total != number_of_keys(mapped_reads)) {cout << "error in classify\n"; exit(1);}
}


void gene::print(){
	for(vector<transcript>::iterator i=begin(); i!=end(); i++)
	{
		cout << i->name << "\t" << i->len << endl;
	}
	cout << "gene has " << mapped_reads.size() << " hits for " << number_of_keys(mapped_reads) << " reads" << endl;
};

inline static int _fill_map( const bam1_t *b, void *d ) {

	bamdata *bd = (bamdata*)d;

	stringstream rname;
	for ( int k = 0; k < b->core.l_qname; k++ ) rname << b->data[k];
	rname << "|" <<  0;

	int pos = b->core.pos;
	int isize = abs(b->core.isize);

//	if(b->core.flag == 16){
		//cout << b->core.flag << endl;
		
	((multimap<string, hit>*)(bd->m))->insert(make_pair(rname.str(),hit(bd->t, pos, isize)));	
//				}
//cout<< rname.str() << endl;
	return 0;
}

void gene::import_reads(vector<bamfile> &bf)
{
	mapped_reads.clear();

	for(int n=0; n<bf.size(); ++n)
	{
		samfile_t *_fs = bf[n].bam;
		bam_index_t *_in = bf[n].bamindex;
		 
		int tid,t0, tlen;

		for ( vector<transcript>::iterator t=begin(); t!=end(); t++ ) {
			bam_parse_region( _fs->header, t->name.c_str(), &tid, &t0, &tlen);
		if (tid != -1){
				bamdata d = {&(*t), &mapped_reads };
				bam_fetch( _fs->x.bam, _in, tid, 0, tlen, &d, _fill_map );
			}
		}
	}
};

#define PERIOD 8
double logprod(vector<double>&x, unsigned int nx)
{
    if(nx == 0) return 0.;
	int exponent=0;
	double prod=1.0;

	vector<double>::iterator it=x.begin();
    prod = frexp(*it,&exponent); 
	++it;

	for(unsigned int n=1; n<nx; ++n, ++it)
	{
		//prod = __frexp(prod*(*it),&expo); 
		if(n % PERIOD ==0) {
			int expo;
			prod = frexp(prod*(*it),&expo); 
			exponent += expo;
		}
		else {
			prod *= (*it); 
		}
	}
	return (log(prod)+exponent)*log(2.0);
}

void sumreads(vector<vector<hit> > &reads, vector<double> &sums, int readlen, vector<double> &frag_distribution)
{
	vector<vector<hit> >::iterator row;
	vector<hit>::iterator col;
	
	int n=0;
	for(row=reads.begin(); row != reads.end(); ++row){
		double sum = 0.;
		for(col = (*row).begin(); col!= (*row).end(); ++col){
			sum += col->trans->prob; //*frag_distribution[col->isize - readlen - 1];
		}
		// ll += log(sum);
	
		sums[n]=sum;
		++n;
	}
}

double gene::compute_likelihood_paired(double *p)
{
	setprobs(p);

	double norm=0.;
	double lnorm = 0.;
	for(vector<transcript>::iterator t=begin(); t!=end(); ++t){norm += t->prob * t->len;}
	lnorm = log(norm);
//	cout << lnorm << endl;
	double ll=0.;

	vector<double> sums(reads.size(),0.);
	sumreads(reads, sums, readlen, frag_distribution);
	
#if(0)
	for(unsigned int n=0; n< sums.size(); ++n)
	{
		// ll += icsi_log(sum,logtable,n_mantissa);
		ll += log(sums[n]);
	}
#endif

#if(1)
	ll=logprod(sums, sums.size());
#endif
//cout << ll << endl;	
	return (ll-reads.size()*lnorm);
};

double gene::compute_likelihood_single(double *p)
//fast version based on number_per_class
//probably doesn't work if we add weights, e.g. positional bias, etc.
{
	setprobs(p);

	 double norm=0.;
	for(vector<transcript>::iterator t=begin(); t!=end(); ++t)
	{
		norm += t->prob * (t->len-readlen);
//		cout << t->weight << endl;
	}
	double lognorm=log(norm);

	// cout << "Norm: " << norm << endl;	
	double ll=0.;	

	for(map<transcript_index, unsigned int>::iterator m_it=number_per_class.begin(); m_it!=number_per_class.end(); ++m_it)
	{
		transcript_index theClass = (*m_it).first;
 
		// cout << (*m_it).first << ":" << (*m_it).second << endl;
		double sum=0.;
		for (unsigned int n=0; n<size(); ++n)
		{
			unsigned int maps_on_n = theClass & 1;
			sum += maps_on_n * at(n).prob;
			theClass = theClass >> 1;
		}
		// cout << sum << endl;
		unsigned int num = (*m_it).second;
		ll += num*(log(sum)-lognorm);
	}
	return ll;
};


#ifndef BAMFILES_HH
#define BAMFILES_HH

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "samtools/sam.h"
#include "tokens.hh"
#include <algorithm>

using namespace std;

typedef struct
{
	char name[512];
	samfile_t *bam;
	bam_index_t *bamindex;
}bamfile;

class bamfiles: public map<int, vector<bamfile> >
{
public:
	bamfiles(char **filesname, unsigned int nf, char* replica)
	{

		vector<string> repnames;
		if(strlen(replica)==0)
		{
			for(unsigned int n=0; n<nf; ++n){
			char temp[512];
			sprintf(temp,"%d",n);
			repnames.push_back(string(temp));
							}
		}
		else	{
		TokenizeV(string(replica),repnames,":");
		}
		
		for(unsigned int n = 0; n<nf; ++n)
		{	
			bamfile bam;
			strcpy(bam.name, filesname[n]);
			bam.bam = samopen(bam.name, "rb", 0 );
			if(!bam.bam){exit(1);}
			bam.bamindex = bam_index_load(bam.name);
			if(!bam.bamindex){exit(1);}
			
			(*this)[atoi(repnames[n].c_str())].push_back(bam);
		}

	};
	
	
};

#endif

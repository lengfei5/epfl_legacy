#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <cstdlib>


#include "tokens.hh"

void Tokenize(const string& str,
	list<string>& tokens,
	const string& delimiters)
{
// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
// Found a token, add it to the list.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void TokenizeV(const string& str,
	vector<string>& tokens,
	const string& delimiters)
{
// Skip delimiters at beginning

	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
// Find first "non-delimiter".

	string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
// Found a token, add it to the list.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void clean_left(string &s)
{	
	while(s.substr(0,1) == " " || s.substr(0,1) == "{"){
		// std::cout << "cl:"  << s << endl;
		s=s.erase(0,1);
	}
		// std::cout << "cl:"  << s << endl;	
}
void erode_left(string &s, char *c)
{	
	if (s.substr(0,1) == c) s=s.erase(0,1);
}

void clean_right(string &s)
{
	while(s.substr(s.size()-1,1) == " " || s.substr(s.size()-1,1) == "}"){
		s=s.erase(s.size()-1,1);
	}
}

void erode_right(string &s, char *c)
{
	if(s.substr(s.size()-1,1) ==c) s=s.erase(s.size()-1,1);
}

void clean_lr(string &s)
{
	clean_right(s);
	clean_left(s);
}

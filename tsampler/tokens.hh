#ifndef TOKENS_HH
#define TOKENS_HH

#include <vector>
#include <list>
#include <string>

using namespace std;

void Tokenize(const string& str,
	list<string>& tokens,
	const string& delimiters);

void TokenizeV(const string& str,
	vector<string>& tokens,
	const string& delimiters);

void clean_left(string &s);
void clean_right(string &s);
void clean_lr(string &s);
void erode_left(string &s, char *c);	// remove one character
void erode_right(string &s, char *c);	//

#endif

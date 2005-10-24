#include "ngila.h"
#include <math.h>
#include <getopt.h>

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

bool g_bNoCase = true;
bool g_bNegate = false;
bool g_bMsg = false;

void print_aln(const string& n1, const string& s1, const string& n2, const string& s2, const string& msg);

double g_dMatch = 0.0;
double g_dReplacement = 1.0;

char g_csUsage[] =
"ngila [hnpisa:b:c:x:m:r:] file\n" \
"  -a -b -c gap parameters\n" \
"  -m match cost\n" \
"  -r mismatch/replacement cost\n" \
"  -x substitution matix\n" \
"  -n negate matrix\n" \
"  -p don't negate matrix\n" \
"  -i case insensitive sequences\n" \
"  -s case sensitive sequence\n" \
"  -q no message\n" \
"  -v message\n" \
"  -h usage\n" \
"\n";

int main(int argc, char *argv[])
{
	int ch;
	char *csMatrix = NULL;
	while((ch = getopt(argc, argv, "hnpisvqa:b:c:x:m:r:")) != -1)
	{
		switch(ch)
		{
		case 'a':
			dA = atof(optarg);
			break;
		case 'b':
			dB = atof(optarg);
			break;
		case 'c':
			dC = atof(optarg);
			break;
		case 'm':
			g_dMatch = atof(optarg);
			break;
		case 'r':
			g_dReplacement = atof(optarg);
			break;
		case 'x':
			csMatrix = optarg;
			break;
		case 'n':
			g_bNegate = true;
			break;
		case 'p':
			g_bNegate = false;
			break;
		case 's':
			g_bNoCase = false;
			break;
		case 'i':
			g_bNoCase = true;
			break;
		case 'q':
			g_bMsg = false;
			break;
		case 'v':
			g_bMsg = true;
			break;
		case 'h':
		case '?':
		default:
			printf("%s", g_csUsage);
			return 1;
			break;
		}
	}
    argc -= optind;
    argv += optind;
	if(csMatrix == NULL)
	{
		for(int i = 0;i<128;++i)
			for(int j=0;j<128;++j)
				mCost[i][j] = ((i == j) ? g_dMatch : g_dReplacement);
	}
	else if(!parse_matrix(csMatrix))
	{
		fprintf(stderr, "Parsing substitution cost matrix failed.\n");
		return 1;
	}
	
	SeqVec seqs;
	StringVec names;
	if(!parse_file(argv[0], names, seqs))
	{
		fprintf(stderr, "Parsing \"%s\" failed.\n", argv[0]);
		return 1;
	}
	for_each(seqs.begin(), seqs.end(), seqproc);
	
	Sequence seqA, seqB;
	ostringstream msg;
	for(size_t i = 0; i < seqs.size(); ++i)
	{
		for(size_t j=i+1; j < seqs.size(); ++j)
		{
			double d = align_pair(seqs[i], seqs[j], seqA, seqB);
			msg.str().clear();
			msg << "Score = " << d;
			print_aln(names[i], seqA, names[j], seqB, msg.str() );
		}
	}


	return 0;
}

char letproc(char ch)
{
	return (g_bNoCase ? (char)toupper(ch) : ch);
}

char letproc(const char *cs)
{
	return letproc(cs[0]);
}

double numproc(const char *cs)
{
	double d = atof(cs);
	return (g_bNegate ? -d : d);
}

void seqproc(string& ss)
{
	// Remove Initial Gap
	string::size_type sy;
	if(ss[0] == '-')
	{
		for(sy = 1; sy < ss.size() && ss[sy] == '-' ; ++sy) { }
		ss.erase(0, sy);
	}
	// Main Body
	for(string::size_type sz = 0; sz < ss.size(); ++sz)
	{
		ss[sz] = letproc(ss[sz]);
		// Remove Gaps
		if(ss[sz+1] == '-')
		{
			for( sy = sz+2; sy < ss.size() && ss[sy] == '-' ; ++sy)
				{ }
			ss.erase(sz+1, sy-sz-1);
		}
	}
}

void print_aln(const string& n1, const string& s1, const string& n2, const string& s2, const string& msg)
{
	cout << "CLUSTAL multiple sequence alignment (Created by Ngila Beta";
	if(g_bMsg)
		cout << ": " << msg;
	cout << ")" << endl << endl << endl;

	size_t sz = s1.size();
	size_t l;
	// Print interleaved sequences
	for(size_t u = 0; u < sz; u+=l)
	{
		l = min(60u, sz);
		// Print a row of each sequence
		cout << setw(15) << setiosflags(ios::left) << n1 << " " << s1.substr(u, l) << endl;
		cout << setw(15) << setiosflags(ios::left) << n2 << " " << s2.substr(u, l) << endl;
		cout << endl << endl;
	}
}

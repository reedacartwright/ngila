%{
#include <stdio.h>
#include <string.h>

#include "ngila.h"

#define strdup _strdup

using namespace std;

int yylex();
void yyerror (char *s);

SeqDB db;

SeqDB::Pos szSeqNum = 0;
SeqDB::Pos szSeq = 0;

string ws2ss(vector<char*>::const_iterator itB, vector<char*>::const_iterator itE );

void delete_words(vector<char*> *vs);

%}

%name-prefix="seq_"
%yacc
%defines
%error-verbose

%start file

%union {
	char ch;
	char *cs;
	std::vector<char*> *vs;
	int  n;
	double d;
}

%token <ch> CLUSTAL
%token <ch>	FASTA
%token <ch> PHYLIP
%token <ch> NEXUS

%token <n>	NUMBER
%token <cs> WORD
%token <ch> ENDL
%token <ch> GT '>'
%token <ch> UNKNOWN

%type <vs>	line
%type <vs>	words
%type <vs>	fsabody
%type <vs>	fsahead

%%
file:
  CLUSTAL alnfile
| FASTA	  fsafile
| PHYLIP  phyfile
| NEXUS	  nexfile
;

alnfile: line alnbody
;

alnbody: alnline 
| alnbody alnline
;

alnline: line {
	if(!$1->empty())
		db.Add(ws2ss($1->begin(), $1->begin()+1), ws2ss($1->begin()+1, $1->end()));
	delete_words($1);
} 
;

fsafile: fsaseq
| fsafile fsaseq
;

fsaseq: fsahead fsabody {
	db.Add(ws2ss($1->begin(), $1->end()), ws2ss($2->begin(), $2->end()));
	delete_words($1);
	delete_words($2);
}
;

fsahead: '>' line { $$ = $2; }
;

fsabody: line { $$ = $1; }
| fsabody line { $$ = $1; $$->insert($$->end(), $2->begin(), $2->end()); delete $2; }
;

phyfile: phyhead phybody
;

phyhead: NUMBER NUMBER ENDL {
	szSeqNum = $1;
	szSeq = 0;
}
;

phybody: phyline
| phybody phyline
;

phyline: line {
	if(!$1->empty())
	{
		if(szSeq < szSeqNum)
			db.Add(ws2ss($1->begin(), $1->begin()+1), ws2ss($1->begin()+1, $1->end()));
		else
			db.Add(szSeq%szSeqNum, ws2ss($1->begin(), $1->end()));
		++szSeq;
	}
	delete_words($1);	
}
;

nexfile:
;

line: words ENDL { $$ = $1; }
| ENDL { $$ = new vector<char*>(); }
;

words: WORD { $$ = new vector<char*>(1, $1); }
| words WORD { $$ = $1; $$->push_back($2); }
;

%%
extern FILE *seq_in;

bool parse_file(const char* csFile, StringVec &vNames, SeqVec &vSeqs)
{
	seq_in = fopen(csFile, "rt");
	if(seq_in == NULL)
		return false;
	db.Clear();
	if(yyparse()!=0)
	{
		db.Clear();
		return false;
	}
	db.Transfer(vNames, vSeqs);	
	return true;
}

string ws2ss(vector<char*>::const_iterator itB, vector<char*>::const_iterator itE )
{
	string ss;
	for(vector<char*>::const_iterator it = itB; it!=itE; ++it)
		ss.append(*it);
	return ss;
}

void delete_words(vector<char*> *vs)
{
	for_each(vs->begin(), vs->end(), delete_func);
	delete vs;
}
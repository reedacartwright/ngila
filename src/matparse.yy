/*  Nigla - Logarithmic Sequence Alignments
    Copyright (C) 2005  Reed A. Cartwright

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
%{
#include "ngila.h"

#ifdef HAVE_STDIO_H
#	include <stdio.h>
#endif

#include <vector>

using namespace std;

int yylex();
void yyerror (char *s);

struct Line
{
	char ch;
	vector<double>* pvd;
};

extern int mat_lineno;

%}

%name-prefix="mat_"
%yacc
%defines
%error-verbose

%start matrix

%union {
	double d;
	char ch;
	Line ln;
	std::vector<double> *pvd;
	std::vector<char> *pvch;
	std::vector<Line> *pvln;
}

%token <d>  	NUM
%token <ch> 	LET
%token <ch> 	UNKNOWN
%token			NL
%token      	END

%type <ln>		line
%type <pvd>		list
%type <pvln>	body
%type <pvch>	header
%type <pvch>	matrix

%%
matrix: header body {
	vector<char> &head = *$1;
	vector<Line> &body = *$2;
	const int n = head.size();
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			mCost[(size_t)body[i].ch][(size_t)head[j]] = (*body[i].pvd)[j];
		}
		delete body[i].pvd;
	}
	delete $1;
	delete $2;
}
;

header: LET { $$ = new vector<char>(1,$1); }
| header LET { $$ = $1; $$->push_back($2); }
;

body: line { $$ = new vector<Line>(1, $1); }
| body line { $$ = $1; $$->push_back($2); }
;

line: LET list { $$.ch = $1; $$.pvd = $2; }
;

list: NUM { $$ = new vector<double>(1, $1); } 
| list NUM { $$ = $1; $$->push_back($2); }
;

%%

extern FILE *mat_in;
extern int mat_lineno;
extern char mat_text[];

bool parse_matrix(const char* csFile)
{
	mat_in = fopen(csFile, "rt");
	if(mat_in == NULL)
		return false;
	return (yyparse()==0);
}

void mat_error(char *s)
{
	fprintf(stderr, "ALERT (line %d): %s: \"%s\".\n", mat_lineno, s, mat_text);
}
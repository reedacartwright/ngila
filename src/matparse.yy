%{
// cmdline: bison -y -d  -o "$(InputName).cc" "$(InputPath)"
//			move  /Y "$(InputName).hh" "$(InputName).h"
// outputs: "$(InputName).cc"; "$(InputName).h"

#include <stdio.h>
#include <vector>

#include "ngila.h"

using namespace std;

int yylex();
void yyerror (char *s);

struct Line
{
	char ch;
	vector<double>* pvd;
};

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

%token <d>  NUM
%token <ch> LET
%token <ch> UNKNOWN
%token		NL
%token      END

%type <ln> line
%type <pvd> list
%type <pvln> body
%type <pvch> header

%%
matrix: header body {
	vector<char> &head = *$1;
	vector<Line> &body = *$2;
	const int n = head.size();
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			mCost[body[i].ch][head[j]] = (*body[i].pvd)[j];
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
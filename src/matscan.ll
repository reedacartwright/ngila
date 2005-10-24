%{
#include <vector>

struct Line
{
	char ch;
	std::vector<double>* pvd;
};

#include "ngila.h"

#include "matparse.h"

using namespace std;

#define yylval mat_lval

int mat_lineno = 1;

%}

%option nounput
%option noyywrap
%option prefix="mat_" outfile="lex.yy.c"

DIGIT  [0-9]
NUMBER [-+]?{DIGIT}+("."{DIGIT}+)?([eE][+-]?{DIGIT}+)?
SPACE [ \t\r\v\f]

%%

"#"[^\n]+ {
	// Skip comments
}

{NUMBER} {
	yylval.d = numproc(yytext);
	return NUM;
}

[A-Za-z*] {
	yylval.ch = letproc(yytext);
	return LET;
}

{SPACE}+ {
	// Skip space
}

\n { mat_lineno++;}

<<EOF>> {
	yyterminate();
}

. {
	yylval.ch = yytext[0];
	return UNKNOWN;	
}

%%

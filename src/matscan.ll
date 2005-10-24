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

%}

%option nounput
%option noyywrap
%option prefix="mat_" outfile="lex.yy.c"
%option yylineno

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

\n { }

<<EOF>> {
	yyterminate();
}

. {
	yylval.ch = yytext[0];
	return UNKNOWN;	
}

%%

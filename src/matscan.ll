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

#include <vector>

struct Line
{
	char ch;
	std::vector<double>* pvd;
};

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

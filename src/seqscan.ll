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
	#include "seqparse.h"
	
	#define yylval seq_lval
	using namespace std;
	int nLine = 1;
%}

%option noyywrap
%option nounput
%option prefix="seq_" outfile="lex.yy.c"

DIGIT	[0-9]
NUMBER	{DIGIT}+
SPACE	[ \t\r\v\f]
NONSP	[^ \t\r\v\f\n>]
WORD	{NONSP}+

%x aln
%x fsa
%x phy
%x nex

%%

"CLUSTAL" {
	yyless(0);
	BEGIN(aln);
	return CLUSTAL;
}

">" {
	yyless(0);
	BEGIN(fsa);
	return FASTA;
}

{SPACE}*{DIGIT}+{SPACE}+{DIGIT}+ {
	yyless(0);
	BEGIN(phy);
	return PHYLIP;
}

"#NEXUS" {
	yyless(0);
	BEGIN(nex);
	return NEXUS;
}

<fsa,aln,phy,nex><<EOF>> {
	BEGIN(INITIAL);
	return ENDL;
}

<fsa>">" {
	yylval.ch = yytext[0];
	return yytext[0];
}

<phy,nex>{NUMBER} {
	yylval.n = atoi(yytext);
	return NUMBER;
}

<*>{WORD} {
	yylval.cs = strdup(yytext);
	return WORD;
}

<*>{SPACE}+ {
	// Skip space
}

<*>\n {
	nLine++;
	return ENDL;
}

<<EOF>> {
	yyterminate();	
}

%%

void seq_error(char *s)
{
	fprintf(stderr, "ALERT (line %d): %s: \"%s\".\n", nLine, s, yytext);
}

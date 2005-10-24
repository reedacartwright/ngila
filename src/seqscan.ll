%{
	#include "ngila.h"
	#include "seqparse.h"
	
	#define yylval seq_lval
	using namespace std;
	int nLine = 1;
%}

%option noyywrap
%option nounput
%option prefix="seq_"

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

. {
	yylval.ch = yytext[0];
	return UNKNOWN;	
}

%%

void seq_error(char *s)
{
	fprintf(stderr, "ALERT (line %d): %s: \"%s\".\n", nLine, s, yytext);
}

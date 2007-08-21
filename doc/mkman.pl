#!/usr/bin/perl -w

use strict;

my $cmds = '';
while(<>)
{
	if(/^\s*XCMD/)
	{
		my $m = $_;
		$m =~ s/^\s*XCMD\(//;
		my @m = split(/\s*,\s*/, $m);
		$m[0] =~ s![)][(]!-!g;
		$m[0] =~ s![()]!!g;
		$m[1] =~ s![()]!!g;
		$m[2] =~ s![\"\']!!g;
		my $arg = ($m[3] =~ /bool/) ? '' : ' arg';
		my $sm = ($m[1] =~ /\w/) ? " [\\-$m[1]]" : '';
		$cmds .= ".TP\n.B \\-\\-$m[0]$sm$arg\n$m[2]\n";
	}
}

print <<EOF;
.TH "Ngila" 1
.SH NAME
Ngila \- Global Pairwise Alignments with Logarithmic and Affine Gap Costs
.SH EXAMPLE
ngila -m zeta -t 0.1 -k 2.0 -r 0.05 -z 1.65 sequences.fasta
.SH DESCRIPTION
Ngila is a global, pairwise sequence alignment program that implements the algorithm of
Miller and Myers (1988), optimized for gap logarithmic and affine gap costs: C(x)=a+b*x+c*ln x.  These costs are more biologically realistic and more powerful than standard affine gap costs.  In addition Ngila implements two generalized pair-HMM models of indel formation, which transform evolutionary parameters to alignment costs.
.SH REFERENCE
Cartwright RA (2007) Ngila: global pairwise alignments with logarithmic and affine gap costs. Bioinformatics. 23(11):1427\-1428
.SH WEBSITE
http://scit.us/projects/ngila/
.SH OPTIONS
$cmds
.SH INPUT
Input sequences must be in FASTA format.
EOF


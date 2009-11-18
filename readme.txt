NGILA VERSION 1.2 - Logarithmic and Affine Sequence Alignments

Copyright (C) 2005-2007  Reed A. Cartwright - All rights reserved.

DESCRIPTION
  Ngila is a global alignment program that can align pairs of sequences using
  logarithmic and affine gap penalties.

REFERENCE
  Cartwright RA (2007) Ngila: global pairwise alignments with logarithmic and
  affine gap costs. Bioinformatics. 23(11):1427-1428

CONTACT
  See AUTHORS for contact information.

LICENSE
  See COPYING for license information.

INSTALLATION
  See INSTALL for installation instructions.

DOWNLOAD
  Ngila can be downloaded from the url <http://scit.us/projects/ngila/>, which
  is its development website.

COMMAND LINE USAGE
  ngila -m zeta -t 0.1 -k 2.0 -r 0.05 -z 1.65 sequences.fasta

  See 'man ngila' or 'ngila --help' for complete command line usage.

  See <http://scit.us/projects/ngila/> for more details on running ngila.

INPUT FILES
  The input file has to be in FASTA format.  If more than two sequences are
  given then Ngila will align based on the 'pairs' option.

OUTPUT FILES
  Output can be in either Clustal or Fasta format and to either the stdout or
  a file.  Clustal is the default.  The format is read from the output file's
  extension or specified directly; "ngila -o seqs.fas" and "ngila -o
  fas:seqs.txt" both produce fasta output.  "ngila -o aln:-" sends Clustal
  formated sequence to stdout.
    

SUBSTITUTION MATRIX
  An example of the format can be seen in matrix/dna

ALGORITHM
  Ngila implements a Miller and Myers (1988) candidate list method of sequence
  alignment with the gap cost being of the form g(x) = a + b*x + c*ln x.  Ngila
  will return the alignment with the minimum cost and has rules for breaking
  ties.  Ngila's main alignment algorithm is divide-and-conquer, which requires
  O(M) memory; but slower than a holistic, O(MN) memory algorithm.

  Ngila implements a secondary, holistic algorithm for alignment, which is
  faster.  The options -M and -N (-M is for the larger sequence) allow users to
  specify thresholds for when the holistic algorithm is used instead of the DnC
  algorithm.  For example, command 'ngila -M 5000 -N 5000 seqs.aln' will align
  the sequences in 'seq.aln' via the divide-and-conquer algorithm, but when
  subsequences less than or equal to 5000-5000 are being aligned, the holistic
  algorithm will be used.


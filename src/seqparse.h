/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     CLUSTAL = 258,
     FASTA = 259,
     PHYLIP = 260,
     NEXUS = 261,
     NUMBER = 262,
     WORD = 263,
     HWORD = 264,
     ENDL = 265,
     GT = 266,
     UNKNOWN = 267
   };
#endif
/* Tokens.  */
#define CLUSTAL 258
#define FASTA 259
#define PHYLIP 260
#define NEXUS 261
#define NUMBER 262
#define WORD 263
#define HWORD 264
#define ENDL 265
#define GT 266
#define UNKNOWN 267




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 59 "seqparse.yy"
typedef union YYSTYPE {
	char ch;
	char *cs;
	std::vector<char*> *vs;
	int  n;
	double d;
} YYSTYPE;
/* Line 1403 of yacc.c.  */
#line 70 "/usr/home/reed/Projects/ngila/current/src/seqparse.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE seq_lval;




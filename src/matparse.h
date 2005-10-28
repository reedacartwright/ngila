/* A Bison parser, made from /usr/home/reed/ngila/current/src/matparse.yy, by GNU bison 1.75.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

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
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

#ifndef BISON_MATPARSE_H
# define BISON_MATPARSE_H

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     LET = 259,
     UNKNOWN = 260,
     NL = 261,
     END = 262
   };
#endif
#define NUM 258
#define LET 259
#define UNKNOWN 260
#define NL 261
#define END 262




#ifndef YYSTYPE
#line 49 "matparse.yy"
typedef union {
	double d;
	char ch;
	Line ln;
	std::vector<double> *pvd;
	std::vector<char> *pvch;
	std::vector<Line> *pvln;
} yystype;
/* Line 1281 of /usr/local/share/bison/yacc.c.  */
#line 63 "matparse.h"
# define YYSTYPE yystype
#endif

extern YYSTYPE mat_lval;


#endif /* not BISON_MATPARSE_H */


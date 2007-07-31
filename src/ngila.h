/****************************************************************************
 *  Copyright (C) 2005-2007  Reed A. Cartwright, PhD <reed@scit.us>         *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                         * 
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

#ifndef NGILA_H
#define NGILA_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

//#define BOOST_SPIRIT_DEBUG 1

#include <iostream>
#include <vector>
#include <map>
#include <string>

// Utility function
#define CERROR(err) ((std::cerr << "ERROR: " << err << endl), false);

// Typedefs
//typedef std::string Sequence;
//typedef std::vector<Sequence> SeqVec;
//typedef std::vector<std::string> StringVec;
//typedef std::map<std::string, Sequence> SeqMap;
//
//// External variables
//extern double mCost[128][128]; //substitution scoring
//extern double dB,dA,dC,dF,dG,dH;
//extern bool g_bFreeEnds;
//extern size_t g_szM, g_szN;
//
//// Static Variables
//static const char chGap = '-';
//
//// Functions
//double align_pair(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD); // align A and B, placing results in C and D
//double align_pair_x(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD);
//char letproc(const char *cs);
//double numproc(const char *cs);
//void seqproc(std::string& ss);
//
//bool parse_matrix(const char* csFile);
//bool parse_file(const char* csFile, StringVec &vNames, SeqVec &vSeqs);
//
//// Templates Functions
//template<typename val> val min3(val a, val b, val c)
//{
//	return std::min(a,std::min(b,c));
//}
//template<typename val, typename rep> val min3(val a, val b, val c, rep& r)
//{
//	if(a <= c)
//		if(a <= b)
//			{ r = (rep)0; return a; }
//		else
//			{ r = (rep)1; return b; }
//	else if(b <= c)
//		{ r = (rep)1; return b; }
//	else
//		{ r = (rep)2; return c; }
//}
//



#endif //NGILA_H

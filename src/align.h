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
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

#ifndef ALIGN_H
#define ALIGN_H

#include <deque>
#include <iostream>
#include <iomanip>
#include <ios>

#include "models.h"
#include "seqdb.h"

const char chGap = '-';

class aligner;

class alignment
{
public:
	typedef seq_db::sequence sequence;
	typedef seq_db::sequence name;
	typedef seq_db::value_type pair_type;
	typedef int aln_atom;
	typedef std::deque<aln_atom> aln_data;
	
	alignment() { }
	alignment(const pair_type &A, const pair_type &B) { set_seqs(A,B); }
	alignment(const sequence &A, const sequence &B) { set_seqs(A,B); }
	
	inline void set_names(const name &A, const name &B) { nameA = A; nameB = B; }
	inline void set_seqs(const sequence &A, const sequence &B) { seqA = A; seqB = B;}
	inline void set_seqs(const pair_type &A, const pair_type &B) { set_names(A.first, B.first); set_seqs(A.second, B.second); }
	
	template<class OS>
	inline void print(OS &os, const char *msg=NULL) const;
	
protected:
	sequence seqA, seqB;
	name nameA, nameB;
	
	aln_data data;
	
	friend class aligner;
};

class aligner
{
public:
	struct indel
	{
		indel() { }
		indel(size_t pp, size_t xx, double dd)
			: p(pp), x(xx), d(dd) { }
		size_t p; // Column or Row
		size_t x; // Crossing Point
		double d; // Score
	};	
	
	typedef alignment::sequence sequence;
	typedef alignment::aln_data aln_data;
	typedef std::vector<indel> indel_vec;
	typedef std::vector<int> travel_row;
	typedef std::vector< travel_row > travel_table;
	
	aligner(const cost_model& m, size_t ma, size_t mb, bool fe) : costs(m), szMa(ma), szMb(mb), bFreeEnds(fe) { };
	
	double align(alignment &rAln);
	
protected:
	cost_model costs;
	size_t szMa;
	size_t szMb;
	bool bFreeEnds;
	
	travel_table tabTravel;
	
private:
	inline double indel_cost(const indel& in, size_t q) const { return in.d + GC[q-in.p]; }
	inline double indel_cost_x(const indel& in) const { return in.d + GC[in.x-in.p]; }
	inline double indel_fcost(const indel& in, size_t q) const { return in.d + FGC[q-in.p]; }
	inline double indel_fcost_x(const indel& in) const { return in.d + FGC[in.x-in.p]; }
	inline double indel_rcost(const indel& in, size_t q) const { return in.d + GC[in.p-q]; }
	inline double indel_rcost_x(const indel& in) const { return in.d + GC[in.p-in.x]; }

	double align_x(const sequence &seqa, const sequence &seqB, aln_data &rAln);
	
	double align_mn(sequence::const_iterator itA1, sequence::const_iterator itA2,
					sequence::const_iterator itB1, sequence::const_iterator itB2,
					aln_data &rAln, bool bFreeFront, bool bFreeBack);
	
	void update_ins_forward(indel_vec &T, size_t i, size_t j, size_t szZ);
	void update_del_forward(indel_vec &T, size_t i, size_t j, size_t szZ);
	void update_del_forward_f(indel_vec &T, size_t i, size_t j, size_t szZ);
	
	std::vector<double> CC[2];
	std::vector<double> RR[2];
	std::vector<double> GC, FGC;
	
	indel_vec T;
	std::vector<indel_vec> SF;
	std::vector<indel_vec> SR;	
};

template<class OS>
inline void alignment::print(OS &os, const char *msg) const
{
	std::string strA, strB, strC;
	sequence::const_iterator ait = seqA.begin(), bit = seqB.begin();
	for(aln_data::const_iterator cit = data.begin(); cit != data.end(); ++cit)
	{
		if(*cit == 0)
		{
			strC.append(1, *ait==*bit ? '*' : ' ');
			strA.append(ait, ait+1); ++ait;
			strB.append(bit, bit+1); ++bit;
		}
		else if(*cit > 0)
		{
			strC.append(*cit, ' ');
			strA.append(*cit, chGap);
			strB.append(bit, bit+*cit); bit+=*cit;
		}
		else
		{
			strC.append(-*cit, ' ');
			strA.append(ait, ait-*cit); ait-=*cit;
			strB.append(-*cit, chGap);
		}
	}
			
	os << "CLUSTAL multiple sequence alignment (Created by " << PACKAGE_STRING;
	if(msg != NULL)
		os << "; " << msg;
	os << ")" << std::endl << std::endl;

	size_t sz = strA.size();
	// Print interleaved sequences
	size_t a=0, b=0;
	std::string ss;
	for(size_t u = 0; u < sz; u+=60)
	{
		// Print a row of each sequence
		ss = strA.substr(u, 60);
		a += ss.length() - std::count(ss.begin(), ss.end(), chGap);
		os << std::setw(14) << std::setiosflags(std::ios::left) << nameA << " " << ss << " " << a << std::endl;
		
		ss = strB.substr(u, 60);
		b += ss.length() - std::count(ss.begin(), ss.end(), chGap);
		os << std::setw(14) << std::setiosflags(std::ios::left) << nameB << " " << ss << " " << b << std::endl;
		
		os << std::setw(14) << std::setiosflags(std::ios::left) << " "   << " " << strC.substr(u, 60) << std::endl;
		os << std::endl;
	}
}

#endif

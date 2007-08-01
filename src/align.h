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

#include "models.h"
#include "seqdb.h"

const char chGap = '-';

class alignment
{
public:
	typedef seq_db::sequence sequence;
	typedef int aln_atom;
	typedef std::deque<aln_atom> aln_data;
	
	inline void swap_seqs()
	{
		swap(seqA, seqB);
		for(aln_data::iterator it = data.begin(); it != data.end(); ++it)
			*it = -*it;
	}
	inline void clear()
	{
		data.clear();
	}
	inline void set_seqs(const sequence &A, const sequence &B)
	{
		seqA = A;
		seqB = B;
	}
	inline void push_front(aln_atom a)
	{
		data.push_front(a);
	}
	inline void push_back(aln_atom a)
	{
		data.push_back(a);
	}
	
	inline void print()
	{
		std::cout << data << std::endl;
	}
	
protected:
	sequence seqA;
	sequence seqB;
	
	aln_data data;
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
	typedef std::vector<indel> indel_vec;
	typedef std::vector<int> travel_row;
	typedef std::vector< travel_row > travel_table;
	
	aligner(const cost_model& m, size_t ma, size_t mb, bool fe) : costs(m), szMa(ma), szMb(mb), bFreeEnds(fe) { };
	
	double align(const sequence &seqA, const sequence &seqB, alignment &rAln);
	
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
		
	double align_x(const sequence &seqA, const sequence &seqB, alignment &rAln);
	double align_mn(sequence::const_iterator itA1, sequence::const_iterator itA2,
					sequence::const_iterator itB1, sequence::const_iterator itB2,
					alignment &rAln, bool bFreeFront, bool bFreeBack);
	
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

#endif

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

#ifdef _MSC_VER
#	pragma warning(disable: 4288)
#endif

#include <algorithm>

#include "ngila.h"
#include "align.h"

using namespace std;

#ifdef min
#	undef min
#endif

template<class T>
inline T min3(T a, T b, T c)
{
	return min(a, min(b,c));
}

inline bool dle(double a, double b)
{
	return (a <= b*(1.0+DBL_EPSILON));
}

double aligner::align(alignment &aln)
{
	double d;
	aln.data.clear();	
	if(aln.seqA.dna.size() >= aln.seqB.dna.size())
		d = align_x(aln.seqA.dna, aln.seqB.dna, aln.data);
	else
	{
		d = align_x(aln.seqB.dna, aln.seqA.dna, aln.data);
		for(aln_data::iterator it = aln.data.begin(); it != aln.data.end(); ++it)
			*it = -*it;
	}
	return d;
}

double aligner::align_x(const sequence &seqA, const sequence &seqB, aln_data &rAln)
{
	size_t sz = seqB.size()+1;
	CC[0].resize(sz);
	CC[1].resize(sz);
	RR[0].resize(sz);
	RR[1].resize(sz);
	SF.resize(sz);
	SR.resize(sz);
	DM.resize(sz);
	GC.resize(1u+max(seqA.size(),seqB.size()));
	FGC.resize(GC.size());
	for(size_t u = 0;u<GC.size();++u)
	{
		GC[u] = costs.gapcost(u);
		FGC[u] = costs.freegapcost(u);
		
	}
	
	// Allocate travel table to maximum possible size
	tabTravel.resize(std::min(szMa,seqA.size())+1);
	size_t szB = std::min(szMb,seqB.size())+1;
	for(travel_table::iterator it = tabTravel.begin();
	    it != tabTravel.end(); ++it)
		it->resize(szB);
	
	//run recursive algorithm

	return align_r(seqA.begin(), seqA.end(), seqB.begin(), seqB.end(), rAln, bFreeEnds, bFreeEnds);
}

double aligner::align_mn(sequence::const_iterator itA1, sequence::const_iterator itA2,
				 sequence::const_iterator itB1, sequence::const_iterator itB2,
				 aln_data &rAln, bool bFreeFront, bool bFreeBack)
{
	// O(MN) memory algorithm
	size_t szNa = itA2-itA1;
	size_t szNb = itB2-itB1;
	
	CC[0][0] = 0.0;
	tabTravel[0][0] = 0;
	for(size_t j=1;j<=szNb;++j)
	{
		CC[0][j] = GC[j];
		tabTravel[0][j] = static_cast<int>(j);
	}
	// TODO: Is this needed?
	SF[0].assign(1, indel(0, szNa, CC[0][0]));

	for(size_t i=1;i<=szNa;++i)
	{
		CC[1][0] = (bFreeFront ? FGC : GC)[i];
		tabTravel[i][0] = -static_cast<int>(i);
		for(size_t j=1;j<=szNb;++j)
		{
			double dM = CC[0][j-1]+costs.mCost[(size_t)itA1[i-1]][(size_t)itB1[j-1]];
			update_ins_forward(T,i,j,szNb);
			double dI = indel_cost(T.back(),j);
			double dD;
			if(bFreeBack && j == szNb)
			{
				update_del_forward_f(SF[j],i,j,szNa);
				dD = indel_fcost(SF[j].back(),i);
			}
			else
			{
				update_del_forward(SF[j],i,j,szNa);
				dD = indel_cost(SF[j].back(),i);
			}
				
			if(dM < dI && dM < dD)
			{
				CC[1][j] = dM;
				tabTravel[i][j] = 0;
			}
			else if(dI <= dM && dI < dD)
			{
				CC[1][j] = dI;
				tabTravel[i][j] = static_cast<int>(j-T.back().p);
			}
			else
			{
				CC[1][j] = dD;
				tabTravel[i][j] = -static_cast<int>(i-SF[j].back().p);
			}
		}
		swap(CC[0], CC[1]);
	}
	size_t i = szNa, j = szNb;

	while(!(i == 0 && j == 0))
	{
		int t = tabTravel[i][j];
		rAln.push_front(t);
		if(t == 0)
		{
			--i;
			--j;
		}
		else if(t > 0)
			j -= static_cast<size_t>(t);
		else
			i -= static_cast<size_t>(-t);
	}
	return CC[0][szNb];
}

double aligner::align_s(sequence::const_iterator itA1, sequence::const_iterator itA2,
				 sequence::const_iterator itB1, sequence::const_iterator itB2,
				 aln_data &rAln, bool bFreeFront, bool bFreeBack)
{
	size_t szNa = itA2-itA1;
	size_t szNb = itB2-itB1;
	if(szNa == 0 && szNb == 0)
		return 0.0; // Nothing to do
	else if(szNb == 0)
	{
		rAln.push_back(-static_cast<alignment::aln_atom>(szNa));
		return (bFreeFront||bFreeBack ? FGC : GC)[szNa]; // delete A
	}
	else if(szNa == 0)
	{
		rAln.push_back(static_cast<alignment::aln_atom>(szNb));
		return GC[szNb]; // insert B
	}
	else if(szNa == 1 && szNb == 1)
	{
		double d1 = costs.mCost[(size_t)itA1[0]][(size_t)itB1[0]];
		double d2 = (bFreeFront||bFreeBack) ? FGC[1]+GC[1] : GC[1]+GC[1];
		if(d1 <= d2)
		{
			rAln.push_back(0);
			return d1;
		}
		else if(bFreeFront)
		{
			rAln.push_back(-1);
			rAln.push_back(1);
		}
		else
		{
			rAln.push_back(1);
			rAln.push_back(-1);
		}
		return d2;
	}
	else if(szNa == 1)
	{
		double dTemp = costs.mCost[(size_t)itA1[0]][(size_t)itB1[0]] + GC[szNb-1];
		double dMin = dTemp;
		size_t i = 0;
		for(size_t j=1;j<szNb-1;++j)
		{
			dTemp = GC[j] + costs.mCost[(size_t)itA1[0]][(size_t)itB1[j]] + GC[szNb-j-1];
			if(dTemp < dMin)
			{
				dMin = dTemp;
				i = j;
			}
		}
		dTemp = costs.mCost[(size_t)itA1[0]][(size_t)itB1[szNb-1]] + GC[szNb-1];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			i = szNb-1;
		}
		dTemp = (bFreeFront ? FGC : GC)[1] + GC[szNb];
		if(dTemp < dTemp)
		{
			dMin = dTemp;
			i = static_cast<size_t>(-1);
		}

		dTemp = (bFreeBack ? FGC : GC)[1] + GC[szNb];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			i = szNb;
		}

		if(i == static_cast<size_t>(-1))
		{
			// Del A, Ins B(1,Nb)
			rAln.push_back(-1);
			rAln.push_back(static_cast<alignment::aln_atom>(szNb));
		}
		else if( i == 0)
		{
			// A->B(1), Ins B(2,Nb)
			rAln.push_back(0);
			rAln.push_back(static_cast<alignment::aln_atom>(szNb-1));
		}		
		else if(i == szNb-1)
		{
			// Ins B(1,N-1), A->B(Nb)
			rAln.push_back(static_cast<alignment::aln_atom>(szNb-1));
			rAln.push_back(0);
		}
		else if( i == szNb)
		{
			// Ins B(1,N), Del A
			rAln.push_back(static_cast<alignment::aln_atom>(szNb));
			rAln.push_back(-1);
		}
		else
		{
			//Ins B(1,i), A->B(i+1), Ins B(i+2,Nb)
			rAln.push_back(static_cast<alignment::aln_atom>(i));
			rAln.push_back(0);
			rAln.push_back(static_cast<alignment::aln_atom>(szNb-i-1));
		}
		return dMin;
	}
	else if(szNb == 1)
	{

		double dTemp = costs.mCost[(size_t)itA1[0]][(size_t)itB1[0]]
				+ (bFreeBack ? FGC : GC)[szNa-1];
		double dMin = dTemp;
		size_t i = 0;
		for(size_t j=1;j<szNa-1;++j)
		{
			dTemp = costs.mCost[(size_t)itA1[j]][(size_t)itB1[0]]
					+ (bFreeFront ? FGC : GC)[j]
					+ (bFreeBack ? FGC : GC)[szNa-j-1];
			if(dTemp < dMin)
			{
				dMin = dTemp;
				i = j;
			}
		}
		dTemp = costs.mCost[(size_t)itA1[szNa-1]][(size_t)itB1[0]]
				+ (bFreeFront ? FGC : GC)[szNa-1];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			i = szNa-1;
		}

		dTemp =  GC[1] + (bFreeBack ? FGC : GC)[szNa];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			i = (size_t)-1;
		}
		dTemp = (bFreeFront ? FGC : GC)[szNa] + GC[1];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			i = szNa;
		}
		
		if(i == (size_t)-1)
		{
			// Ins B, Del A(1,Na)
			rAln.push_back(1);
			rAln.push_back(-static_cast<alignment::aln_atom>(szNa));
		}
		else if( i == 0)
		{
			// B->A(1), Del A(2,Na)
			rAln.push_back(0);
			rAln.push_back(-static_cast<alignment::aln_atom>(szNa-1));
		}		
		else if(i == szNa-1)
		{
			// Del A(1,Na-1), B->A(Na)
			rAln.push_back(-static_cast<alignment::aln_atom>(szNa-1));
			rAln.push_back(0);
		}
		else if( i == szNa)
		{
			// Del A(1,Na), Ins B
			rAln.push_back(-static_cast<alignment::aln_atom>(szNa));
			rAln.push_back(1);
		}
		else
		{
			//Del A(1,i), B->A(i+1), Del A(i+2,Na)
			rAln.push_back(-static_cast<alignment::aln_atom>(i));
			rAln.push_back(0);
			rAln.push_back(-static_cast<alignment::aln_atom>(szNa-i-1));
		}
		return dMin;
	}
	return 0.0;
}


double aligner::align_r(sequence::const_iterator itA1, sequence::const_iterator itA2,
				 sequence::const_iterator itB1, sequence::const_iterator itB2,
				 aln_data &rAln, bool bFreeFront, bool bFreeBack)
{
	size_t szNa = itA2-itA1;
	size_t szNb = itB2-itB1;
	size_t szMh = szNa/2;

	if(szNa <= 1 || szNb <= 1) // Handle Special Cases
		return align_s(itA1, itA2, itB1, itB2, rAln, bFreeFront, bFreeBack);
	else if(szNa <= szMa && szNb < szMb)
		return align_mn(itA1, itA2, itB1, itB2, rAln, bFreeFront, bFreeBack);

	CC[0][0] = 0.0;
	for(size_t j=1;j<=szNb;++j)
		CC[0][j] = GC[j];
	
	SF[0].assign(1, indel(0, szNa, CC[0][0]));
	// Forward Algorithm
	for(size_t i=1;i<=szMh;++i)
	{
		CC[1][0] = (bFreeFront ? FGC : GC)[i];
		for(size_t j=1;j<=szNb;++j)
		{
			update_ins_forward(T,i,j,szNb);
			update_del_forward(SF[j],i,j,szNa);
			CC[1][j] = min3(CC[0][j-1]+costs.mCost[(size_t)itA1[i-1]][(size_t)itB1[j-1]],
				indel_cost(T.back(),j), indel_cost(SF[j].back(),i));
		}
		swap(CC[0], CC[1]);
	}
	//CC[0] is now C(M/2,j)

	//Reverse Algorithm
	RR[0][szNb] = 0.0;
	DM[szNb].c = bFreeBack ? indel_fcost(SF[szNb][0],szNa)
					: indel_cost(SF[szNb][0],szNa);
	DM[szNb].s = 0;
	DM[szNb].z = 0;
	DM[szNb].x = szNa;
	RR[0][0] = GC[szNb];
	DM[0].c = indel_cost(SF[0][0],szNa)+RR[0][0];
	DM[0].s = 0;
	DM[0].z = 0;
	DM[0].x = szNa;
	for(size_t j=szNb-1;j>0;--j)
	{
		RR[0][j] = GC[szNb-j];
		DM[j].c = indel_cost(SF[j][0],szNa)+RR[0][j];
		DM[j].s = 0;
		DM[j].z = 0;
		DM[j].x = szNa;
	}
	for(size_t i=szNa-1;i!=szMh-1;--i)
	{
		if( SF[szNb].size() < DM[szNb].z+1 && i <= SF[szNb][DM[szNb].z+1].x)
			++DM[szNb].z; // Advance position
		RR[1][szNb] = (bFreeBack ? FGC : GC)[szNa-i];
		double dTemp = indel_cost(SF[szNb][DM[szNb].z],i)+RR[1][szNb];
		if(dTemp < DM[szNb].c)
		{
			DM[szNb].c = dTemp;
			DM[szNb].s = DM[szNb].z;
			DM[szNb].x = i;
		}
			
		for(size_t j=szNb-1;j!=(size_t)-1;--j)
		{
			update_ins_reverse(T,i,j,szNb);
			update_del_reverse(SR[j],i,j,szNa);
			// Type I
			double d1 = RR[0][j+1]+costs.mCost[(size_t)itA1[i]][(size_t)itB1[j]];
			double d2 = indel_rcost(T.back(),j);
			double d3 = indel_rcost(SR[j].back(),i);
			RR[1][j] = min3(d1,d2,d3);

			// Minimum Type II cost
			if( SF[j].size() < DM[j].z+1 && i <= SF[j][DM[j].z+1].x)
				++DM[j].z; // Advance position
			if(bFreeFront && j==0)
				dTemp = indel_fcost(SF[j][DM[j].z],i)+RR[1][j];
			else
				dTemp = indel_cost(SF[j][DM[j].z],i)+RR[1][j];
			if(dTemp < DM[j].c)
			{
				DM[j].c = dTemp;
				DM[j].s = DM[j].z;
				DM[j].x = i;
			}
		}
		swap(RR[0], RR[1]);
	}

	double dMin = CC[0][0]+RR[0][0];
	size_t pp = szMh;
	size_t xx = pp;
	size_t jj = 0;
	if(dle(DM[0].c, dMin))
	{
		dMin = DM[0].c;
		pp = SF[0][DM[0].s].p;
		xx = DM[0].x;
	}
	for(size_t j=1;j<=szNb;++j)
	{
		double dTemp = CC[0][j]+RR[0][j];
		if(dTemp < dMin)
		{
			dMin = dTemp;
			pp = szMh;
			xx = pp;
			jj = j;
		}
		dTemp = DM[j].c;
		if(dle(dTemp,dMin))
		{
			dMin = dTemp;
			pp = SF[j][DM[j].s].p;
			xx = DM[j].x;
			jj = j;
		}
	}	
	align_r(itA1, itA1+pp, itB1, itB1+jj, rAln, bFreeFront, false);
	if(xx != pp)
		rAln.push_back(-static_cast<alignment::aln_atom>(xx-pp)); // Delete itA1+pp .. itA1+xx
	align_r(itA1+xx, itA2, itB1+jj, itB2, rAln, false, bFreeBack);
	return dMin;	
}

void aligner::update_ins_forward(indel_vec &T, size_t /*i*/, size_t j, size_t szZ)
{
	if(j == 1)
	{
		T.clear();
		T.push_back(indel(0, szZ, CC[1][0]));
		return;
	}
	if( j > T.back().x)
		T.pop_back();
	if( CC[1][j-1] + GC[1] < indel_cost(T.back(), j) ||
		CC[1][j-1] + GC[szZ-j+1] <= indel_cost(T.back(), szZ) )
	{
		while(T.size() && (CC[1][j-1] + GC[T.back().x - j+1] < indel_cost_x(T.back()) ||
			  CC[1][j-1] + GC[szZ-j+1] <= indel_cost(T.back(), szZ)))
			T.pop_back();
		T.push_back(indel(j-1, (T.size() ?
			(T.back().p+costs.kstar(j-1-T.back().p, CC[1][j-1]-T.back().d))
			: szZ) ,CC[1][j-1]));
	}
}

void aligner::update_del_forward(indel_vec &T, size_t i, size_t j, size_t szZ)
{
	if(i == 1)
	{
		T.clear();
		T.push_back(indel(0, szZ, CC[0][j]));
		return;
	}
	if( i > T.back().x)
		T.pop_back();
	if( CC[0][j] + GC[1] < indel_cost(T.back(), i) ||
		CC[0][j] + GC[szZ-i+1] <= indel_cost(T.back(),szZ) )
	{
		while(T.size() && (CC[0][j] + GC[T.back().x - i+1] < indel_cost_x(T.back()) ||
			  CC[0][j] + GC[szZ-i+1] <= indel_cost(T.back(), szZ)))
			T.pop_back();
		T.push_back(indel(i-1, (T.size() ?
			(T.back().p+costs.kstar(i-1-T.back().p, CC[0][j]-T.back().d))
			: szZ), CC[0][j]));
	}
}

void aligner::update_ins_reverse(indel_vec &T, size_t i, size_t j, size_t szZ)
{
	if(j == szZ-1)
	{
		T.clear();
		T.push_back(indel(szZ, 0, RR[1][szZ]));
		return;
	}
	if( j < T.back().x)
		T.pop_back();
	if( RR[1][j+1] + GC[1] <  indel_rcost(T.back(),j) ||
		RR[1][j+1] + GC[j+1] <= indel_rcost(T.back(),0) )
	{
		while(T.size() && (RR[1][j+1] + GC[j+1-T.back().x] < indel_rcost_x(T.back()) ||
			  RR[1][j+1] + GC[j+1] <= indel_rcost(T.back(),0)))
			T.pop_back();
		T.push_back(indel(j+1, (T.size() ?
			(T.back().p-costs.kstar(T.back().p-(j+1), RR[1][j+1]-T.back().d))
			: 0) ,RR[1][j+1]));
	}
}

void aligner::update_del_forward_f(indel_vec &T, size_t i, size_t j, size_t szZ)
{
	if(i == 1)
	{
		T.clear();
		T.push_back(indel(0, szZ, CC[0][j]));
		return;
	}
	if( i > T.back().x)
		T.pop_back();
	if( CC[0][j] + FGC[1] < indel_fcost(T.back(), i) ||
		CC[0][j] + FGC[szZ-i+1] <= indel_fcost(T.back(),szZ) )
	{
		while(T.size() && (CC[0][j] + FGC[T.back().x - i+1] < indel_fcost_x(T.back()) ||
			  CC[0][j] + FGC[szZ-i+1] <= indel_fcost(T.back(), szZ)))
			T.pop_back();
		T.push_back(indel(i-1, (T.size() ?
			(T.back().p+costs.kstar_f(i-1-T.back().p, CC[0][j]-T.back().d))
			: szZ), CC[0][j]));
	}
}



void aligner::update_del_reverse(indel_vec &T, size_t i, size_t j, size_t szZ)
{
	if(i == szZ-1)
	{
		T.clear();
		T.push_back(indel(szZ, 0, RR[0][j]));
		return;
	}
	if( i < T.back().x)
		T.pop_back();
	if( RR[0][j] + GC[1] <  indel_rcost(T.back(),i) ||
		RR[0][j] + GC[i+1] <= indel_rcost(T.back(),0) )
	{
		while(T.size() && (RR[0][j] + GC[i+1-T.back().x] < indel_rcost_x(T.back()) ||
			  RR[0][j] + GC[i+1] <= indel_rcost(T.back(),0)))
			T.pop_back();
		T.push_back(indel(i+1, (T.size() ?
			(T.back().p-costs.kstar(T.back().p-(i+1), RR[0][j]-T.back().d))
			: 0) ,RR[0][j]));
	}
}


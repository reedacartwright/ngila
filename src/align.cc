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

#include "ngila.h"

#include <algorithm>

#include <float.h>
#include <math.h>

using namespace std;

double mCost[128][128];

double dB = 1.0, dA = 1.0, dC = 0.0;

struct MinMidCost
{
	double c; // Minimum cost
	size_t s; // minimum position in s-vector
	size_t x; // bottom of minimum indel
	size_t z; // active position in s-vector
};

vector<MinMidCost> DM;

inline double gapcost(size_t L)
{
	return (L != 0) ? dA+dB*L+dC*log((double)L) : 0.0;
}

inline double kstar(double x, double y)
{
	return x/(1.0-exp((-y+dB*x)/dC)); //needs to handle dB == 0.0?
}

vector<double> CC[2];
vector<double> RR[2];
vector<double> GC;

class Indel
{
public:
	Indel() { }
	Indel(size_t pp, size_t xx, double dd)
		: p(pp), x(xx), d(dd) { }
	size_t p; // Column or Row
	size_t x; // Crossing Point
	double d; // Score

	double Cost(size_t q) { return d + GC[q-p]; }
	double CostX() { return d + GC[x-p]; }
	double RCost(size_t q) { return d + GC[p-q]; }
	double RCostX() { return d + GC[p-x]; }
};

typedef vector<Indel> IndelVec;

IndelVec T;
vector<IndelVec> SF;
vector<IndelVec> SR;

void update_ins_forward(IndelVec& T, size_t i, size_t j, size_t szZ);
void update_del_forward(IndelVec& T, size_t i, size_t j, size_t szZ);
void update_ins_reverse(IndelVec& T, size_t i, size_t j, size_t szZ);
void update_del_reverse(IndelVec& T, size_t i, size_t j, size_t szZ);
double align_pair_r(Sequence::const_iterator itA1, Sequence::const_iterator itA2,
				 Sequence::const_iterator itB1, Sequence::const_iterator itB2,
				 Sequence& seqC, Sequence& seqD,bool bFreeFront, bool bFreeBack);

inline size_t g1(size_t p, size_t x, size_t j, size_t m, size_t n)
{
	size_t r;
	if(p != j)
		r = min(p,j);
	else if(p != x)
		r = p;
	else if(n != m )
		r = min(n,m);
	else
		r = n;
	if(m-x != n-j)
		r += max(p,j)+x-p+max(m-x,n-j);
	else if(p != x)
		r += max(p,j)+x-p;
	else if( p != j)
		r += max(p,j);
	else
		r += n;
	return r;
}

inline size_t g2(size_t p, size_t x, size_t j, size_t m, size_t n)
{
	size_t r;
	if(j <= p)
		r = x+x;
	else
		r = j+j+x-p;
	if(m-x != n-j)
		r += min(m-x,n-j);
	
	return 0;
}

double align_pair(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD)
{
	if(seqA.size() >= seqB.size())
		return align_pair_x(seqA, seqB, seqC, seqD);
	else
		return align_pair_x(seqB, seqA, seqD, seqC);
}

double align_pair_x(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD)
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
	for(size_t u = 0;u<GC.size();++u)
		GC[u] = gapcost(u);
	//run recursive algorithm
	seqC.clear();
	seqD.clear();

	return align_pair_r(seqA.begin(), seqA.end(), seqB.begin(), seqB.end(), seqC, seqD, g_bFreeEnds, g_bFreeEnds);
}

double align_pair_r(Sequence::const_iterator itA1, Sequence::const_iterator itA2,
				 Sequence::const_iterator itB1, Sequence::const_iterator itB2,
				 Sequence& seqA, Sequence& seqB, bool bFreeFront, bool bFreeBack)
{
	size_t szM = itA2-itA1;
	size_t szN = itB2-itB1;
	size_t szMh = szM/2;

	if(szN == 0 && szM == 0)
	{
		return 0.0; // Nothing to do
	}
	else if(szN == 0)
	{
		seqA.append(itA1, itA2);
		seqB.append(szM, chGap);
		return (bFreeFront||bFreeBack) ? 0.0 : GC[szM]; // delete A
	}
	else if(szM == 0)
	{
		seqA.append(szN, chGap);
		seqB.append(itB1, itB2);
		return (bFreeFront||bFreeBack) ? 0.0 : GC[szN]; // insert B
	}
	else if(szM == 1 && szN == 1)
	{
		double d1 = mCost[(size_t)itA1[0]][(size_t)itB1[0]];
		double d2 = bFreeFront ? 0.0 : GC[1];
		d2 += bFreeBack ? 0.0 : GC[1];
		if( d1 <= d2)
		{
			seqA.append(1, itA1[0]);
			seqB.append(1, itB1[0]);
			return d1;
		}
		else if(itA1[0] <= itB1[0])
		{
			seqA.append(1, chGap);
			seqA.append(1, itA1[0]);
			seqB.append(1, itB1[0]);
			seqB.append(1, chGap);
		}
		else
		{
			seqA.append(1, itA1[0]);
			seqA.append(1, chGap);
			seqB.append(1, chGap);
			seqB.append(1, itB1[0]);
		}
		return d2;
	}
	
	if(szM == 1)
	{
		double dTemp = mCost[(size_t)itA1[0]][(size_t)itB1[0]]
					+(bFreeBack ? 0.0 : GC[szN-1]);
		double dMin = dTemp;
		size_t i = 0;
		for(size_t j=1;j<szN;++j)
		{
			dTemp = mCost[(size_t)itA1[0]][(size_t)itB1[j]]
					+(bFreeFront ? 0.0 : GC[j])+(bFreeBack ? 0.0 : GC[szN-j-1]);
			if(dTemp < dMin)
			{
				dMin = dTemp;
				i = j;
			}
		}
		dTemp = (bFreeFront ? 0.0 : GC[1])+(bFreeBack ? 0.0 : GC[szN]);
		if(dTemp < dMin)
		{
			i = (size_t)-1;
			dMin = dTemp;
		}
		dTemp = (bFreeFront ? 0.0 : GC[szN])+(bFreeBack ? 0.0 : GC[1]);
		if(dTemp < dMin)
		{
			i = szN;
			dMin = dTemp;
		}
		if(i == (size_t)-1)
		{
			// Del A, Ins B(1,N)
			seqA.append(1, itA1[0]);
			seqA.append(szN, chGap);
			seqB.append(1, chGap);				
		}
		else if( i == 0)
		{
			// A->B(1), Ins B(2,N)
			seqA.append(1, itA1[0]);			
			seqA.append(szN-1, chGap);
		}		
		else if(i == szN-1)
		{
			// Ins B(1,N-1), A->B(N)
			seqA.append(szN-1, chGap);			
			seqA.append(1, itA1[0]);			
		}
		else if( i == szN)
		{
			// Ins B(1,N), Del A
			seqA.append(szN, chGap);
			seqA.append(1, itA1[0]);
			seqB.append(1, chGap);						
		}
		else
		{
			//Ins B(1,i), A->B(i+1), Ins B(i+2,N)
			seqA.append(i, chGap);
			seqA.append(1, itA1[0]);
			seqA.append(szN-i-1, chGap);	
		}
		seqB.append(itB1, itB2);
		return dMin;
	}
	else if(szN == 1)
	{
		double dTemp = mCost[(size_t)itA1[0]][(size_t)itB1[0]]
					+(bFreeBack ? 0.0 : GC[szM-1]);
		double dMin = dTemp;
		size_t i = 0;
		for(size_t j=1;j<szM;++j)
		{
			dTemp = mCost[(size_t)itA1[j]][(size_t)itB1[0]]
					+(bFreeFront ? 0.0 : GC[j])+(bFreeBack ? 0.0 : GC[szM-j-1]);
			if(dTemp < dMin)
			{
				dMin = dTemp;
				i = j;
			}
		}
		dTemp = (bFreeFront ? 0.0 : GC[1])+(bFreeBack ? 0.0 : GC[szM]);
		if(dTemp < dMin)
		{
			i = (size_t)-1;
			dMin = dTemp;
		}
		dTemp = (bFreeFront ? 0.0 : GC[szM])+(bFreeBack ? 0.0 : GC[1]);
		if(dTemp < dMin)
		{
			i = szM;
			dMin = dTemp;
		}
		if(i == (size_t)-1)
		{
			// Del B, Ins A(1,M)
			seqB.append(1, itB1[0]);
			seqB.append(szM, chGap);
			seqA.append(1, chGap);				
		}
		else if( i == 0)
		{
			// B->A(1), Ins A(2,M)
			seqB.append(1, itB1[0]);			
			seqB.append(szM-1, chGap);
		}		
		else if(i == szM-1)
		{
			// Ins A(1,M-1), B->A(M)
			seqB.append(szM-1, chGap);			
			seqB.append(1, itB1[0]);			
		}
		else if( i == szM)
		{
			// Ins A(1,M), Del B
			seqB.append(szM, chGap);
			seqB.append(1, itB1[0]);
			seqA.append(1, chGap);						
		}
		else
		{
			//Ins A(1,i), B->A(i+1), Ins A(i+2,M)
			seqB.append(i, chGap);
			seqB.append(1, itB1[0]);
			seqB.append(szM-i-1, chGap);	
		}
		seqA.append(itA1, itA2);
		return dMin;
	}

	CC[0][0] = 0.0;
	for(size_t j=1;j<=szN;++j)
		CC[0][j] = bFreeFront ? 0.0 : GC[j];
	
	SF[0].push_back(Indel(0, szM, CC[0][0]));
	// Foward Algorithm
	for(size_t i=1;i<=szMh;++i)
	{
		CC[1][0] = bFreeFront ? 0.0 : GC[i];
		for(size_t j=1;j<=szN;++j)
		{
			update_ins_forward(T,i,j,szN);
			update_del_forward(SF[j],i,j,szM);
			CC[1][j] = min3(CC[0][j-1]+mCost[(size_t)itA1[i-1]][(size_t)itB1[j-1]],
				T.back().Cost(j), SF[j].back().Cost(i));
		}
		swap(CC[0], CC[1]);
	}
	//CC[0] is now C(M/2,j)

	//Reverse Algorithm
	RR[0][szN] = 0.0;
	DM[szN].s = 0;
	DM[szN].z = 0;
	DM[szN].x = szM;
	DM[szN].c = bFreeBack ? SF[j][0].d : SF[j][0].Cost(szM);
	for(size_t j=szN;j!=(size_t)-1;--j)
	{
		RR[0][j] = bFreeBack ? 0.0 : GC[szN-j];
		DM[j].s = 0;
		DM[j].z = 0;
		DM[j].x = szM;
		DM[j].c = SF[j][0].Cost(szM)+RR[0][j];
	}
	for(size_t i=szM-1;i!=szMh-1;--i)
	{
		RR[1][szN] = bFreeBack ? 0.0 : GC[szM-i];
		if( SF[j].size() < DM[j].z+1 && i <= SF[j][DM[j].z+1].x)
		  ++DM[j].z; // Advance position
		double dTemp = SF[j][DM[j].z].Cost(i)+RR[1][j];
		if(dTemp < DM[j].c)
			{
				DM[j].c = dTemp;
				DM[j].s = DM[j].z;
				DM[j].x = i;
			}
			
		for(size_t j=szN-1;j!=(size_t)-1;--j)
		{
			update_ins_reverse(T,i,j,szN);
			update_del_reverse(SR[j],i,j,szM);
			// Type I
			double d1 = RR[0][j+1]+mCost[(size_t)itA1[i]][(size_t)itB1[j]];
			double d2 = T.back().RCost(j);
			double d3 = SR[j].back().RCost(i);
			RR[1][j] = min3(d1,d2,d3);
			// Minimum Type II cost
			if( SF[j].size() < DM[j].z+1 && i <= SF[j][DM[j].z+1].x)
				++DM[j].z; // Advance position
			double dTemp = SF[j][DM[j].z].Cost(i)+RR[1][j];
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
	printf("%f: %d %d %d\n", dMin, pp, xx, jj);
	size_t gg1 = g1(pp,xx,jj, szM, szN);
	size_t gg2 = g2(pp,xx,jj, szM, szN);
	size_t g = g1(SF[0][DM[0].s].p,DM[0].x,0, szM, szN);
	size_t h = g2(SF[0][DM[0].s].p,DM[0].x,0, szM, szN);
	if(DM[0].c < dMin || (DM[0].c == dMin &&
			(g > gg1 || (g == gg1 &&
			(h < gg2 || (h == gg2 &&
			DM[0].x > xx))))))
	{
		dMin = DM[0].c;
		pp = SF[0][DM[0].s].p;
		xx = DM[0].x;
		gg1 = g;
		gg2 = h;
	}
	printf("%f: %d %d %d\n", DM[0].c,SF[0][DM[0].s].p , DM[0].x, 0);
	for(size_t j=1;j<=szN;++j)
	{
		double dTemp = CC[0][j]+RR[0][j];
		printf("%f: %d %d %d\n", dTemp, szMh,szMh, j);
		g = g1(szMh,szMh,j, szM, szN);
		h = g2(szMh,szMh,j, szM, szN);
		if(dTemp < dMin || (dTemp == dMin &&
			(g > gg1 || (g == gg1 &&
			(h < gg2 )))))
		{
			dMin = dTemp;
			pp = szMh;
			xx = pp;
			jj = j;
			gg1 = g;
			gg2 = h;
		}
		dTemp = DM[j].c;
		g = g1(SF[j][DM[j].s].p, DM[j].x, j, szM, szN);
		h = g2(SF[j][DM[j].s].p, DM[j].x, j, szM, szN);
		if(dTemp < dMin || (dTemp == dMin &&
			(g > gg1 || (g == gg1 &&
			(h < gg2 || (h == gg2 &&
			DM[j].x > xx))))))
		{
			dMin = dTemp;
			pp = SF[j][DM[j].s].p;
			xx = DM[j].x;
			jj = j;
			gg1 = g;
			gg2 = h;
		}
		printf("%f: %d %d %d\n", DM[j].c,SF[j][DM[j].s].p , DM[j].x, j);
	}	
	
	printf("Min %f: %d %d %d\n", dMin, pp, xx, jj);
	
	align_pair_r(itA1, itA1+pp, itB1, itB1+jj, seqA, seqB, bFreeFront, false);
	if(xx != pp)
	{
		// Delete itA1+pp .. itA1+xx
		seqA.append(itA1+pp, itA1+xx);
		seqB.append(xx-pp, chGap);
	}
	align_pair_r(itA1+xx, itA2, itB1+jj, itB2, seqA, seqB, false, bFreeBack);
	return dMin;	
}

void update_ins_forward(IndelVec& T, size_t /*i*/, size_t j, size_t szZ)
{
	if(j == 1)
	{
		T.clear();
		T.push_back(Indel(0, szZ, CC[1][0]));
		return;
	}
	if( j > T.back().x)
		T.pop_back();
	if( CC[1][j-1] + GC[1] < T.back().Cost(j) ||
		CC[1][j-1] + GC[szZ-j+1] <= T.back().Cost(szZ) )
	{
		while(T.size() && (CC[1][j-1] + GC[T.back().x - j+1] < T.back().CostX() ||
			  CC[1][j-1] + GC[szZ-j+1] <= T.back().Cost(szZ)))
			T.pop_back();
		T.push_back(Indel(j-1, (T.size() ?
			(T.back().p+(size_t)kstar(j-1-T.back().p, CC[1][j-1]-T.back().d))
			: szZ) ,CC[1][j-1]));
	}
}

void update_del_forward(IndelVec& T, size_t i, size_t j, size_t szZ)
{

	if(i == 1)
	{
		T.clear();
		T.push_back(Indel(0, szZ, CC[0][j]));
		return;
	}
	if( i > T.back().x)
		T.pop_back();
	if( CC[0][j] + GC[1] < T.back().Cost(i) ||
		CC[0][j] + GC[szZ-i+1] <= T.back().Cost(szZ) )
	{
		while(T.size() && (CC[0][j] + GC[T.back().x - i+1] < T.back().CostX() ||
			  CC[0][j] + GC[szZ-i+1] <= T.back().Cost(szZ)))
			T.pop_back();
		T.push_back(Indel(i-1, (T.size() ?
			(T.back().p+(size_t)kstar(i-1-T.back().p, CC[0][j]-T.back().d))
			: szZ), CC[0][j]));
	}
}

void update_ins_reverse(IndelVec& T, size_t /*i*/, size_t j, size_t szZ)
{
	//i; // i should be constant in all comparisons
	if(j == szZ-1)
	{
		T.clear();
		T.push_back(Indel(szZ, 0, RR[1][szZ]));
		return;
	}
	if( j < T.back().x)
		T.pop_back();
	if( RR[1][j+1] + GC[1] <  T.back().RCost(j) ||
		RR[1][j+1] + GC[j] <= T.back().RCost(1) )
	{
		while(T.size() && (RR[1][j+1] + GC[j+1-T.back().x] < T.back().RCostX() ||
			  RR[1][j+1] + GC[j] <= T.back().RCost(1)))
			T.pop_back();
		T.push_back(Indel(j+1, (T.size() ?
			(T.back().p-(size_t)kstar(T.back().p-j-1, RR[1][j+1]-T.back().d)+1)
			: 0) ,RR[1][j+1]));
	}
}

void update_del_reverse(IndelVec& T, size_t i, size_t j, size_t szZ)
{
	if(i == szZ-1)
	{
		T.clear();
		T.push_back(Indel(szZ, 0, RR[0][j]));
		return;
	}
	if( i < T.back().x)
		T.pop_back();
	if( RR[0][j] + GC[1] <  T.back().RCost(i) ||
		RR[0][j] + GC[i] <= T.back().RCost(1) )
	{
		while(T.size() && (RR[0][j] + GC[i+1-T.back().x] < T.back().RCostX() ||
			  RR[0][j] + GC[i] <= T.back().RCost(1)))
			T.pop_back();
		T.push_back(Indel(i+1, (T.size() ?
			(T.back().p-(size_t)kstar(T.back().p-i-1, RR[0][j]-T.back().d)+1)
			: 0) ,RR[0][j]));
	}
}


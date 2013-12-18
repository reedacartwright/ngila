/****************************************************************************
 *  Copyright (C) 2005-2010  Reed A. Cartwright, PhD <reed@scit.us>         *
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

#include "ngila.h"
#include "models.h"

#include <boost/math/special_functions/zeta.hpp>

#include <limits>
#include <iomanip>

using namespace std;

inline double log_zeta(double z) {
	return log(boost::math::zeta<double>(z));
}

/****************************************************************************
 *    class cost_model                                                      *
 ****************************************************************************/


bool cost_model::create(const ngila_app::args &rargs)
{
	dA = rargs.cost_intersection;
	dB = rargs.cost_linear;
	dC = rargs.cost_logarithmic;
	dF = rargs.cost_intersection_free;
	dG = rargs.cost_linear_free;
	dH = rargs.cost_logarithmic_free;

	if(rargs.cost_matrix.empty())
	{		
		for(size_t i=0;i<sub_matrix_size;++i)
			for(size_t j=0;j<sub_matrix_size;++j)
				mCost[i][j] = (i == j) ? rargs.cost_match : rargs.cost_mismatch;
	}
	else
	{
		if(!parse_matrix(rargs.cost_matrix.c_str(), mCost, rargs.case_insensitivity))
			return CERRORR("parsing of \'" << rargs.cost_matrix.c_str() << "\' failed.");
	}
	return true;
}

/****************************************************************************
 *    class k2p_model                                                       *
 ****************************************************************************/

enum {nA, nC, nG, nT, nN, nSize};
const int nuc_table[] = {
	/*64*/	-1, nA, nN, nC, nN, -1, -1, nG, nN, -1, -1, nN, -1, nN, nN, -1,
	/*80*/	-1, -1, nN, nN, nT, nT, nN,	nN, -1, nN, -1, -1, -1, -1, -1, -1,
	/*96*/	-1, nA, nN, nC, nN, -1, -1, nG,	nN, -1, -1, nN, -1, nN, nN, -1,
	/*112*/	-1, -1, nN, nN, nT, nT, nN,	nN, -1, nN, -1, -1, -1, -1, -1, -1,
	/*128*/
};
const int pupy[] = {1, 0, 1, 0, -1};

bool k2p_model::create(const ngila_app::args &rargs)
{
	if(rargs.branch_length <= 0.0)
		return CERRORR("branch length must be positive.");
	if(rargs.ratio <= 0.0)
		return CERRORR("Ts/Tv ratio must be positive.");
	if(rargs.avgaln <= 0.0)
		return CERRORR("avgaln must be positive.");
	if(rargs.indel_rate <= 0.0)
		return CERRORR("indel rate must be positive.");
	
	double p_ts = 0.25-0.5*exp(-rargs.branch_length*(2.0*rargs.ratio+1.0)/(rargs.ratio+1.0))
			+ 0.25*exp(-2.0*rargs.branch_length/(rargs.ratio+1.0));
	double p_tv = 0.5-0.5*exp(-2.0*rargs.branch_length/(rargs.ratio+1.0));
	double p_match = 1.0-p_ts-p_tv;
	double l_h = -2.0*rargs.indel_rate*rargs.branch_length
		+log(rargs.avgaln)-log(rargs.avgaln+1.0);
	
	double c_ts = -(log(p_ts)+l_h+log(0.25));
	double c_tv = -(log(p_tv)+l_h+log(0.125));
	double c_m  = -(log(p_match)+l_h+log(0.25));
	
	dNucScale = rargs.no_scale ? 0.0 : 0.5*c_m;

	double sub_costs[nSize][nSize];
	for(size_t i = 0;i<nN;++i)
	{
		sub_costs[i][i] = c_m;
		for(size_t j = 0; j < i; ++j)
		{
			sub_costs[i][j] = sub_costs[j][i] = (pupy[i] == pupy[j]) ? c_ts : c_tv;
		}
		sub_costs[nN][i] = sub_costs[i][nN] = -(l_h+log(0.0625)); // All N's are weighted by 0.25
	}
	sub_costs[nN][nN] = -(l_h+log(0.0625)); // All N's are weighted by 0.25
	
	sub_matrix_clear(mCost);
	for(size_t i = 0;i<64;++i)
	{
		int ni = nuc_table[i];
		for(size_t j=0;j<64;++j)
		{
			int nj = nuc_table[j];
			if(ni == -1 || nj == -1)
				continue;
			mCost[i+64][j+64] = sub_costs[ni][nj] - 2.0*dNucScale;
		}
	}
	dEnd = rargs.no_scale ? log(rargs.avgaln+1.0) : 0.0;
	
	return true;
}

template<class T>
struct isambnuc
{
	bool operator()(T n)
	{
		return (nuc_table[n-64] == nN);
	}
};

double k2p_model::offset(const std::string &seqA, const std::string &seqB) const
{
	isambnuc<std::string::value_type> cpred;
	size_t ncount = count_if(seqA.begin(), seqA.end(), cpred);
	ncount += count_if(seqB.begin(), seqB.end(), cpred);

	return dEnd-static_cast<double>(ncount)*log(4.0);
}

/****************************************************************************
 *    class zeta_model                                                      *
 ****************************************************************************/

bool zeta_model::create(const ngila_app::args &rargs) {
	if(!k2p_model::create(rargs))
		return false;
	if(rargs.indel_slope <= 1.0)
		return CERRORR("indel slope must be greater than 1.0");
	
	dA = -(log(0.5)+log(1.0-exp(-2.0*rargs.indel_rate*rargs.branch_length))
		+ log(rargs.avgaln)-log(rargs.avgaln+1.0)
		- log_zeta(rargs.indel_slope));
	dB = -log(0.25) - dNucScale; // All N's are weighted by 0.25
	dC = rargs.indel_slope;

	dF = dH = 0.0;
	//dG = -log(0.25) - dNucScale;
	dG = rargs.no_scale ? -log(0.25) : 0.0;

	return true;
}

/****************************************************************************
 *    class geo_model                                                       *
 ****************************************************************************/

bool geo_model::create(const ngila_app::args &rargs) {
	if(!k2p_model::create(rargs))
		return false;
	if(rargs.indel_mean <= 1.0)
		return CERRORR("indel mean must be greater than 1.0");
	
	dA = -(log(0.5)+log(1.0-exp(-2.0*rargs.indel_rate*rargs.branch_length))
		+ log(rargs.avgaln)-log(rargs.avgaln+1.0)
		- log(rargs.indel_mean-1));
	dB = -(log(0.25)+log(rargs.indel_mean-1.0)-log(rargs.indel_mean)) - dNucScale; // All N's are weighted by 0.25
	dC = 0.0;

	dF = dH = 0.0;
	//dG = -log(0.25) - dNucScale;
	dG = rargs.no_scale ? -log(0.25) : 0.0;

	return true;
}

/****************************************************************************
 *    class aa_model                                                        *
 ****************************************************************************/

bool aa_model::create(const ngila_app::args &rargs) {
	enum {D=0,E=20,U=40};
	static const char AA[] = "ARNDCQEGHILKMFPSTWYV";
	static const char aa[] = "arndcqeghilkmfpstwyv";
	
	static const double data[440] = {
#		include "lgmod.incl"
	};
	if(rargs.branch_length <= 0.0)
		return CERRORR("branch length must be positive.");
	if(rargs.avgaln <= 0.0)
		return CERRORR("avgaln must be positive.");
	if(rargs.indel_rate <= 0.0)
		return CERRORR("indel rate must be positive.");
	
	// skeleton parameters
	double l_h = -2.0*rargs.indel_rate*rargs.branch_length
		+log(rargs.avgaln)-log(rargs.avgaln+1.0);	
	
	
	double el[20], x[20][20], m[20];
	sub_matrix_clear(mCost);
	
	// cost of stationary amino acids
	std::fill(&vAACost[0], &vAACost[128],
		std::numeric_limits<double>::max()/16.0);
	for(int i=0;i<20;++i)
		vAACost[AA[i]] = vAACost[aa[i]] = -2.0*log(data[D+i]);
		
	// exp the eigenvalues	
	for(int i=0;i<20;++i)
		el[i] = exp(rargs.branch_length*data[E+i]);
	// D*U*el
	for(int i=0;i<20;++i) {
		for(int j=0;j<20;++j) {
			x[i][j] = data[U+j+20*i]*el[j]/data[D+i];
		}
	}
	// x*(U^T*D^-1)
	// matches
	for(int i=0;i<20;++i) {
		double scost = 0.0;
		for(int k=0;k<20;++k)
			scost += x[i][k]*data[U+k+20*i]*data[D+i];
		scost = -(log(scost))-l_h-vAACost[AA[i]];
		m[i] = 0.0;
		if(!rargs.no_scale) {
			m[i] = scost+2*vAACost[AA[i]];
			scost = -2*vAACost[AA[i]];
		}
		mCost[AA[i]][AA[i]] =  mCost[AA[i]][aa[i]]
			= mCost[aa[i]][AA[i]] =  mCost[aa[i]][aa[i]]
			= scost;
	}
	// mismatches
	for(int i=0;i<20;++i) {
		for(int j=i+1;j<20;++j) {
			double scost = 0.0;
			for(int k=0;k<20;++k)
				scost += x[i][k]*data[U+k+20*j]*data[D+j];
			scost = -(log(scost))-l_h-vAACost[AA[j]];
			scost -= (m[i]+m[j])/2;
			
			mCost[AA[i]][AA[j]] =  mCost[AA[j]][AA[i]]
				= mCost[AA[i]][aa[j]] =  mCost[AA[j]][aa[i]]
				= mCost[aa[i]][AA[j]] =  mCost[aa[j]][AA[i]]
				= mCost[aa[i]][aa[j]] =  mCost[aa[j]][aa[i]]
				= scost;
		}
	}
	dEnd = rargs.no_scale ? log(rargs.avgaln+1.0) : 0.0;
	
	return true;	
}

double aa_model::offset(const string &seqA, const string &seqB) const {
	double off = 0.0;
	for(string::size_type k=0;k<seqA.size();++k) {
		off += vAACost[seqA[k]];
	}
	for(string::size_type k=0;k<seqB.size();++k) {
		off += vAACost[seqB[k]];
	}
	return dEnd+off;
}

/****************************************************************************
 *    class aazeta_model                                                    *
 ****************************************************************************/

bool aazeta_model::create(const ngila_app::args &rargs) {
	if(!aa_model::create(rargs))
		return false;
	if(rargs.indel_slope <= 1.0)
		return CERRORR("indel slope must be greater than 1.0");
	
	dA = -(log(0.5)+log(1.0-exp(-2.0*rargs.indel_rate*rargs.branch_length))
		+ log(rargs.avgaln)-log(rargs.avgaln+1.0)
		- log_zeta(rargs.indel_slope));
	dB = 0.0;
	dC = rargs.indel_slope;
	dF = dH = 0.0;
	dG = 0.0;

	return true;
}

/****************************************************************************
 *    class aageo_model                                                     *
 ****************************************************************************/

bool aageo_model::create(const ngila_app::args &rargs) {
	if(!aa_model::create(rargs))
		return false;
	if(rargs.indel_mean <= 1.0)
		return CERRORR("indel mean must be greater than 1.0");
	
	dA = -(log(0.5)+log(1.0-exp(-2.0*rargs.indel_rate*rargs.branch_length))
		+ log(rargs.avgaln)-log(rargs.avgaln+1.0)
		- log(rargs.indel_mean-1.0));
	dB = -(log(rargs.indel_mean-1.0)-log(rargs.indel_mean));
	dC = 0.0;
	dF = dH = 0.0;
	dG = 0.0;

	return true;
}




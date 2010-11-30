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

#include "ngila.h"
#include "models.h"

#include <limits>

#ifdef USE_GSL_ZETA
#include <gsl/gsl_sf_zeta.h>
inline double log_zeta(double z) { return log(gsl_sf_zeta(z)); }
#else
double log_zeta(double z) {
	// y = log((z-1)*zeta(z)) for z = seq(1,10,0.1)
	static double y[91] = { 0.0,
		0.0568007048928621, 0.111824418961511, 0.165162481245235, 0.216899224977457,
		0.267112722170842, 0.315875427900892, 0.363254740373158, 0.409313489900198,
		0.454110367565165, 0.497700302470746, 0.540134794961249, 0.58146221198276,
		0.621728049750889, 0.660975168080968, 0.699244000065426, 0.736572740229038,
		0.772997513832678, 0.80855252961234, 0.843270217918634, 0.877181355951437,
		0.910315181555989, 0.94269949685307, 0.974360762811216, 1.00532418572837,
		1.03561379646993, 1.0652525232069, 1.09426225830845, 1.12266391996661,
		1.15047750906371, 1.17772216173545, 1.20441619803188, 1.23057716703451,
		1.25622188874908, 1.28136649305977, 1.30602645600061, 1.33021663357389,
		1.35395129332193, 1.37724414383803, 1.40010836238451, 1.42255662076912,
		1.44460110961723, 1.46625356116382, 1.48752527067833, 1.50842711662482,
		1.52896957965087, 1.54916276049038, 1.56901639685795, 1.58853987940593,
		1.60774226680912, 1.62663230003676, 1.64521841586656, 1.66350875969082,
		1.68151119766112, 1.69923332821384, 1.71668249301608, 1.73386578736786,
		1.75079007009434, 1.76746197295883, 1.78388790962537, 1.80007408419733,
		1.81602649935684, 1.83175096412782, 1.84725310128395, 1.86253835442144,
		1.87761199471499, 1.89247912737422, 1.90714469781664, 1.92161349757211,
		1.93589016993295, 1.94997921536273, 1.96388499667612, 1.97761174400137,
		1.99116355953612, 2.00454442210694, 2.01775819154196, 2.03080861286571,
		2.0436993203246, 2.05643384125104, 2.06901559977376, 2.08144792038136,
		2.09373403134592, 2.10587706801283, 2.11788007596313, 2.12974601405368,
		2.1414777573409, 2.15307809989286, 2.16454975749482, 2.17589537025256,
		2.18711750509814, 2.19821865820189
	};
	// truncate at z=10
	if(z >= 10.0)
		return y[90]-log(10.0-1.0);
	// if z is out of range, just return 0.0
	if(z <= 1.0)
		return 0.0;
	// find position in table
	size_t p = static_cast<size_t>((z-1)*10.0);
	// linear approximation
	return y[p]+(10.0*z-(10+p))*(y[p+1]-y[p]) - log(z-1.0);
	
}
#endif

using namespace std;

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
	
	for(size_t i = 0;i<64;++i)
	{
		int ni = nuc_table[i];
		for(size_t j=0;j<64;++j)
		{
			int nj = nuc_table[j];
			if(ni == -1 || nj == -1)
				mCost[i+64][j+64] = numeric_limits<double>::quiet_NaN();
			else
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

bool zeta_model::create(const ngila_app::args &rargs)
{
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

	for(double d=1.001;d<11;d+=0.007) {
		cerr << d << "\t" << log_zeta(d) << endl;
	}

	return true;
}

/****************************************************************************
 *    class geo_model                                                       *
 ****************************************************************************/

bool geo_model::create(const ngila_app::args &rargs)
{
	if(!k2p_model::create(rargs))
		return false;
	if(rargs.indel_mean <= 1.0)
		return CERRORR("indel mean must be greater than 1.0");
	
	dA = -(log(0.5)+log(1.0-exp(-2.0*rargs.indel_rate*rargs.branch_length))
		+ log(rargs.avgaln)-log(rargs.avgaln+1.0)
		- log(rargs.indel_mean));
	dB = -(log(0.25)+log(rargs.indel_mean-1.0)-log(rargs.indel_mean)) - dNucScale; // All N's are weighted by 0.25
	dC = 0.0;
	dF = dH = 0.0;
	//dG = -log(0.25) - dNucScale;
	dG = rargs.no_scale ? -log(0.25) : 0.0;

	return true;
}

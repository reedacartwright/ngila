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

#ifndef MODELS_H
#define MODELS_H

#include <cmath>
#include <cfloat>
#include <limits>

#include <boost/cstdint.hpp>
#include <boost/preprocessor/cat.hpp>

#include "ngila_app.h"
#include "matparser.h"
#include "seqdb.h"

#ifndef TABLE_CELL_BITSIZE
#	define TABLE_CELL_BITSIZE 16
#endif

typedef BOOST_PP_CAT(boost::int, BOOST_PP_CAT(TABLE_CELL_BITSIZE, _t)) _travel_cell;

struct cost_model
{
	double dA, dB, dC;
	double dF, dG, dH;
	sub_matrix mCost;
		
	virtual bool create(const ngila_app::args &rargs);
	virtual double offset(const std::string &seqA,
		const std::string &seqB) const
	{
		return 0.0;
	}
	
	inline double gapcost(size_t L)
	{
		if(L < 1)
			return 0.0;
		if(L > (size_t)std::numeric_limits<_travel_cell>::max())
			return std::numeric_limits<double>::max()/16.0;
		return dA+dB*L+dC*log((double)L);
	}
	inline double freegapcost(size_t L)
	{
		if(L < 1)
			return 0.0;
		if(L > (size_t)std::numeric_limits<_travel_cell>::max())
			return std::numeric_limits<double>::max()/16.0;
		return dF+dG*L+dH*log((double)L);
	}
	inline size_t kstar(double x, double y)
	{
		return static_cast<size_t>(floor(x/(1.0-exp((-y+dB*x)/dC)))); //needs to handle dC == 0.0?
	}
	inline size_t kstar_f(double x, double y)
	{
		return static_cast<size_t>(floor(x/(1.0-exp((-y+dG*x)/dH)))); //needs to handle dC == 0.0?
	}		
};

struct k2p_model : public cost_model {
	double dEnd;
	
	virtual bool create(const ngila_app::args &rargs);
	virtual double offset(const std::string &seqA,
		const std::string &seqB) const;
protected:
	double dNucScale;
};

struct zeta_model : public k2p_model {
	virtual bool create(const ngila_app::args &rargs);
};

struct geo_model : public k2p_model {
	virtual bool create(const ngila_app::args &rargs);
};

struct aa_model : public cost_model {
	double dEnd;
	double vAACost[128];
		
	virtual bool create(const ngila_app::args &rargs);
	virtual double offset(const std::string &seqA,
		const std::string &seqB) const;
protected:
	double dAaScale;
};

struct aazeta_model : public aa_model {
	virtual bool create(const ngila_app::args &rargs);
};

struct aageo_model : public aa_model {
	virtual bool create(const ngila_app::args &rargs);
};


#endif


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


/****************************************************************************
 *    class cost_model                                                      *
 ****************************************************************************/


bool cost_model::create(const ngila_app::args &rargs)
{
	dA = rargs.cost_intersection;
	dB = rargs.cost_linear;
	dC = rargs.cost_logarithmic;
	if(rargs.free_end_gaps)
	{
		dF = rargs.cost_intersection_free;
		dG = rargs.cost_linear_free;
		dH = rargs.cost_logarithmic_free;
	}
	else
	{
		dF = dA;
		dG = dB;
		dH = dC;
	}
	if(rargs.cost_matrix.empty())
	{		
		for(size_t i=0;i<sub_matrix_size;++i)
			for(size_t j=0;j<sub_matrix_size;++j)
				mCost[i][j] = (i == j) ? rargs.cost_match : rargs.cost_mismatch;
	}
	else
	{
		if(!parse_matrix(rargs.cost_matrix.c_str(), mCost, rargs.case_insensitivity))
			return CERROR("parsing of \'" << rargs.cost_matrix.c_str() << "\' failed.");
	}
	return true;
}

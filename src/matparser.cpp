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
#include "matparser.h"

using namespace std;

bool parse_matrix(const char *cs, sub_matrix &rsm)
{
	file_iterator<char> file_first(cs);
	if(!file_first)
	{
		cerror() << "unable to open \'" << cs << "\'" << endl;
		return false;			
	}
	file_iterator<char>  file_last = file_first.make_end();
	
	mat_work w;
	mat_grammar my_grammar(w);
	parse_info< file_iterator<char> > info = parse(file_first, file_last, my_grammar);
	if (!info.full)
	{
		cerror() << "unable to parse \'" << cs << "\'" << endl;
		return false;
	}
	if(!w.process(rsm))
	{
		cerror() << "matrix \'" << cs << "\' specified incorrectly" << endl;
		return false;
	}
	return true;
}

bool mat_work::process(sub_matrix &m)
{
	fill(&m[0][0], &m[128][128], 0.0);

	for(size_t r = 0; r < labels.size(); ++r)
	{
		char rch = labels[r];
		for(size_t c = 0; c < labels.size(); ++c)
		{
			char cch = labels[c];
			m[rch][cch] = data[r][c];
		}
	}
	return false;
}


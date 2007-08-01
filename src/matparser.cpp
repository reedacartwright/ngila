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

bool parse_matrix(const char *cs, sub_matrix &rsm, bool bi)
{
	mat_work w;
	mat_grammar my_grammar(w);
	if(strcmp(cs, "-")==0)
	{
		string ss;
		getline(cin, ss, '\004');
		if(ss.empty())
			return CERROR("unable to open stdin.");
		parse_info<string::const_iterator> info = parse(ss.begin(), ss.end(), my_grammar, blank_p);
		if (!info.full)
			return CERROR("unable to parse stdin.");
	}
	else
	{
		file_iterator<char> file_first(cs);
		if(!file_first)
			return CERROR("unable to open \'" << cs << "\'.");
		file_iterator<char>  file_last = file_first.make_end();
		parse_info< file_iterator<char> > info = parse(file_first, file_last, my_grammar, blank_p);
		if (!info.full)
			return CERROR("unable to parse \'" << cs << "\'.");
	}
	if(!w.process(rsm, bi))
		return CERROR("matrix \'" << cs << "\' must be triangular or square.");
	return true;
}

bool mat_work::process(matrix &m, bool bi) const
{
	sub_matrix_clear(m);
	
	size_t sz = labels.size();

	if(data.size() != sz || !(data.front().size() == sz
		|| data.back().size() == sz) )
		return false;
	labels_type llabels(labels);
	if(bi)
	{
		for(labels_type::iterator it = llabels.begin();
			it != llabels.end(); ++it)
			*it = toupper(*it);
	}
	if(data.front().size() == 1)
	{
		// Assume lower triangular matrix
		for(size_t i = 0; i < sz; ++i)
		{
			if(data[i].size() != i+1)
				return false;
			size_t ich = static_cast<size_t>(llabels[i]);
			for(size_t j = 0; j <= i; ++j)
			{
				size_t jch = static_cast<size_t>(llabels[j]);
				m[ich][jch] = m[jch][ich] = data[i][j];
			}
		}
	}
	else if(data.back().size() == 1)
	{
		// Assume upper triangular matrix
		for(size_t i = 0; i < sz; ++i)
		{
			if(data[i].size() != sz-i)
				return false;
			size_t ich = static_cast<size_t>(llabels[i]);
			for(size_t j = 0; j < sz-i; ++j)
			{
				size_t jch = static_cast<size_t>(llabels[j+i]);
				m[ich][jch] = m[jch][ich] = data[i][j];
			}
		}
	}
	else
	{
		// Assume square matrix
		for(size_t i = 0; i < sz; ++i)
		{
			if(data[i].size() != sz)
				return false;
			size_t ich = static_cast<size_t>(llabels[i]);
			for(size_t j = 0; j < sz; ++j)
			{
				size_t jch = static_cast<size_t>(llabels[j]);
				m[ich][jch] = data[i][j];
			}
		}
	}
	return true;
}



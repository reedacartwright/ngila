/***************************************************************************
 *  Copyright (C) 2005-2007  Reed A. Cartwright, PhD <reed@scit.us>        *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "seqparser.h"

using namespace std;

bool parse_file(const char* cs, seq_db &rdb)
{
	file_iterator<char> file_first(cs);
	if(!file_first)
	{
		cerr << "ERROR: unable to open \'" << cs << "\'" << endl;
		return false;			
	}
	file_iterator<char>  file_last = file_first.make_end();
		
	stack<string> my_stack;
	seq_grammar my_grammar(my_stack, rdb);
	parse_info< file_iterator<char> > info = parse(file_first, file_last, my_grammar);
	if (!info.full)
	{
		cerr << "ERROR: unable to parse \'" << cs << "\'" << endl;
		return false;
	}
	return true;
}


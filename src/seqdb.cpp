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

#include <cstring>
#include <string>
#include <algorithm>
#include <boost/bind.hpp>

#include "seqdb.h"
#include "seqparser.h"

using namespace std;
using namespace boost;

bool seq_db::parse_file(const char *csfile, bool bappend, bool bi, int dirs)
{
	if(!bappend)
		clear();

	stack<string> my_stack;
	seq_grammar my_grammar(my_stack, *this);
	
	if(strcmp(csfile, "-")==0) {
		string ss; 
		getline(cin, ss, '\004'); // EOT, end of transmission, ^D
		if(ss.empty())
			return CERRORR("unable to open stdin.");

		parse_info<string::const_iterator> info = parse(ss.begin(), ss.end(), my_grammar);
		if (!info.full)
			return CERRORR("unable to parse stdin.");
	} else {
		file_iterator<char> file_first(csfile);
		if(!file_first)
			return CERRORR("unable to open \'" << csfile << "\'");
		file_iterator<char>  file_last = file_first.make_end();
		
		parse_info< file_iterator<char> > info = parse(file_first, file_last, my_grammar);
		if (!info.full)
			return CERRORR("unable to parse \'" << csfile << "\'");
	}
		
	for(container::iterator it = cont.begin(); it != cont.end(); ++it)
		cont.modify(it, bind(&seq_data::finalize, _1, bi, dirs));
	
	return true;
}

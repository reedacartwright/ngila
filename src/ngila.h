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
 *                                                                         * 
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ****************************************************************************/

#ifndef NGILA_H
#define NGILA_H

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

//#define BOOST_SPIRIT_DEBUG 1

#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <ostream>
#include <iterator>

#define VERSION_MSG PACKAGE_STRING "\n" \
	"    Copyright (C) 2005-2009  Reed A. Cartwright, PhD <reed@scit.us>\n" \
	"    Report bugs to " PACKAGE_BUGREPORT

// Utility functions
#define CERRORR(err_msg) ((std::cerr << "ERROR: " << err_msg << std::endl), false);
#define CERROR(err_msg) (std::cerr << "ERROR: " << err_msg << std::endl);

namespace std {
template<typename _Tp, typename _CharT, typename _Traits>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& os, const std::vector<_Tp> &v)
{
	if(v.size() == 1)
		os << v.front();
	else if(v.size() > 1)
	{
		std::copy(v.begin(), v.end()-1, std::ostream_iterator<_Tp, _CharT, _Traits>(os, " "));
		os << v.back();
	} 
	return os;
}
template<typename _Tp, typename _CharT, typename _Traits>
basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& os, const std::deque<_Tp> &v)
{
	if(v.size() == 1)
		os << v.front();
	else if(v.size() > 1)
	{
		std::copy(v.begin(), v.end()-1, std::ostream_iterator<_Tp, _CharT, _Traits>(os, " "));
		os << v.back();
	} 
	return os;
}
}

#endif //NGILA_H


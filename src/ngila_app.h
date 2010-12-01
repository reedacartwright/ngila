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

#ifndef NGILA_APP_H
#define NGILA_APP_H
#pragma once

#include "ngila.h"

#include <vector>
#include <utility>

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

/****************************************************************************
 *    class ngila_app                                                       *
 ****************************************************************************/

class ngila_app  {
public:
	typedef std::vector< std::pair<size_t,size_t> > pair_vec;
	
	ngila_app(int argc, char *argv[]);
	virtual ~ngila_app() { }
	virtual int run();

	std::string runname;
	po::options_description desc, indesc;
	po::positional_options_description pdesc;
	po::variables_map vm;

	struct args {			
		// use X-Macros to specify argument variables
		#define XCMD(lname, sname, desc, type, def) type _V(lname) ;
		#include "ngila.cmds"
		#undef XCMD
		std::vector< std::string > input;
	};
	enum { CONST_ALIGN_ORDER = 8, CONST_ALIGN_SEQS = 16 }; 
	
protected:
	args arg;
	
};

namespace boost { namespace program_options {
template<>
inline typed_value<bool>* value(bool* v) {
	return bool_switch(v);
}
template<>
inline typed_value<boost::tribool>* value(boost::tribool* v) {
	typed_value<boost::tribool>* r = new typed_value<boost::tribool>(v);
    r->implicit_value(true, "on");
	return r;
}

// modified from boost/libs/program_options/src/value_semantic.cpp
inline void validate(boost::any& v, const std::vector<std::string>& xs, boost::tribool*, int) {
    check_first_occurrence(v);
	std::string s(get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = char(tolower(s[i]));

    if (s.empty() || s == "on" || s == "yes" || s == "1" || s == "true")
		v = any(boost::tribool(true));
    else if (s == "off" || s == "no" || s == "0" || s == "false")
		v = any(boost::tribool(false));
    else if (s == "null" || s == "maybe" || s == "2" || s == "indeterminate")
		v = any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(invalid_option_value(s));
}
#if !defined(BOOST_NO_STD_WSTRING)
inline void validate(any& v, const std::vector<std::wstring>& xs, boost::tribool*, int)
{
    check_first_occurrence(v);
	std::wstring s(get_single_string(xs, true));

    for (size_t i = 0; i < s.size(); ++i)
        s[i] = wchar_t(tolower(s[i]));

    if (s.empty() || s == L"on" || s == L"yes" || s == L"1" || s == L"true")
        v = any(boost::tribool(true));
    else if (s == L"off" || s == L"no" || s == L"0" || s == L"false")
        v = any(boost::tribool(false));
    else if (s == L"null" || s == L"maybe" || s == L"2" || s == L"indeterminate")
		v = any(boost::tribool(boost::indeterminate));
    else
        boost::throw_exception(invalid_option_value(s));
}
#endif
}}

#endif

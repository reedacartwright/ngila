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

#ifndef NGILA_APP_H
#define NGILA_APP_H

#include "ngila.h"

#include <vector>
#include <utility>

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
	
	po::options_description desc;	
	
	struct args
	{			
	// use X-Macros to specify argument variables
#	define XCMD(lname, sname, desc, type, def) type _V(lname) ;
#	include "ngila.cmds"
#	undef XCMD
	};
	enum { CONST_ALIGN_ORDER = 8, CONST_ALIGN_SEQS = 16 }; 
	
protected:
	args arg;
	
};

#endif


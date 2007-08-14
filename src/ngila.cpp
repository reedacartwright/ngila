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
#include <boost/preprocessor.hpp>
#include <sstream>

#include "ngila_app.h"
#include "seqdb.h"
#include "matparser.h"
#include "models.h"
#include "align.h"

using namespace std;

int main(int argc, char *argv[])
{
	int ret = EXIT_FAILURE;
	try {
		ngila_app app(argc, argv);
		ret = app.run();
	} catch(exception &e) {
		CERROR(e.what());
	}
	return ret;
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value() { return bool_switch(); }

template<>
typed_value<bool>* value(bool* v) { return bool_switch(v); }

}}

ngila_app::ngila_app(int argc, char* argv[]) : desc("Allowed Options")
{
	try {
		desc.add_options()
			#define XCMD(lname, sname, desc, type, def) ( \
				_SS("-", lname) BOOST_PP_EXPR_IF(_SD((_)sname), "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg._JS(_, lname))->default_value(def), \
				desc )
			#include "ngila.cmds"
			#undef XCMD
			;
		po::variables_map vm;
		po::positional_options_description pdesc;
		pdesc.add("input", -1);
		po::store(po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
		po::notify(vm);
		if(!arg.arg_file.empty())
		{
			if(arg.arg_file == "-")
			{
				po::store(po::parse_config_file(cin, desc), vm);	
			}
			else
			{
				std::ifstream ifs(arg.arg_file.c_str());
				if(!ifs.is_open())
				{
					string sse = "unable to open argument file ";
					sse += arg.arg_file;
					throw std::runtime_error(sse);
				}
				po::store(po::parse_config_file(ifs, desc), vm);
			}
			po::notify(vm);
		}
	} catch (exception &e) {
		CERROR(e.what());
		throw std::runtime_error("unable to process command line");
	}
}

int ngila_app::run()
{
	if(arg.version)
	{
		cerr << endl << VERSION_MSG << endl << endl;
		return EXIT_SUCCESS;
	}
	if(arg.help)
	{
		cerr << endl << VERSION_MSG << endl << endl;
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	if(arg.quiet)
	{
		cerr.clear(ios::failbit);
	}
	
	seq_db mydb(arg.remove_gaps);
	for(vector<string>::const_iterator cit = arg.input.begin(); cit != arg.input.end(); ++cit)
	{
		if(!mydb.parse_file(cit->c_str(), true))
		{
			CERROR("parsing of \'" << cit->c_str() << "\' failed.");
			return EXIT_FAILURE;
		}
	}
	cost_model *pmod = NULL;
	if(string("zeta").find(arg.model) == 0)
		pmod = new zeta_model;
	else if(string("geo").find(arg.model) == 0)
		pmod = new geo_model;
	else if(string("cost").find(arg.model) == 0)
		pmod = new cost_model;
	else
	{
		CERROR("unknown model \'" << arg.model << "\'.");
		return EXIT_FAILURE;
	}

	if(!pmod->create(arg))
	{
		CERROR("creating model \'" << arg.model << "\'.");
		return EXIT_FAILURE;
	}
	aligner alner(*pmod, arg.threshold_larger, arg.threshold_smaller, arg.free_end_gaps);
	if(mydb.size() < 2)
	{
		CERROR("two or more sequences are required for alignment.");
		return EXIT_FAILURE;
	}
	for(size_t i=1;i<mydb.size();++i)
	{
		for(size_t j=0;j<i;++j)
		{
			alignment aln(mydb[j], mydb[i]);
			double dcost = alner.align(aln);
			dcost += pmod->offset(mydb[j].second, mydb[i].second);
			ostringstream msg;
			msg << "Cost = " << dcost;
			
			aln.print(cout, msg.str().c_str());
			if(j != mydb.size()-2)
				cout << "//" << endl;
		}
	}
	
	return EXIT_SUCCESS;
}

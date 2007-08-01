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

	seq_db mydb;
	for(vector<string>::const_iterator cit = arg.input.begin(); cit != arg.input.end(); ++cit)
	{
		if(!mydb.parse_file(cit->c_str(), true))
		{
			CERROR("parsing of \'" << cit->c_str() << "\' failed.");
			return EXIT_FAILURE;
		}
	}
	cost_model *pmod = new cost_model;
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
	alignment aln;
	alner.align(mydb[0].second, mydb[1].second, aln);
	aln.print();
	
	return EXIT_SUCCESS;
}



//void print_aln(const string& n1, const string& s1, const string& n2, const string& s2, const char * msg)
//{
//	cout << "CLUSTAL multiple sequence alignment (Created by " << PACKAGE_STRING;
//	if(g_bMsg)
//		cout << ": " << msg;
//	cout << ")" << endl << endl << endl;
//
//	size_t sz = s1.size();
//	size_t l;
//	// Print interleaved sequences
//	for(size_t u = 0; u < sz; u+=l)
//	{
//		l = std::min((size_t)60u, sz);
//		// Print a row of each sequence
//		cout << setw(15) << setiosflags(ios::left) << n1 << " " << s1.substr(u, l) << endl;
//		cout << setw(15) << setiosflags(ios::left) << n2 << " " << s2.substr(u, l) << endl;
//		cout << endl << endl;
//	}
//}

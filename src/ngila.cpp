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

#include <fstream>
#include <ostream>
#include <iterator>

#include "ngila_app.h"
#include "seqdb.h"
#include "matparser.h"

using namespace std;

int main(int argc, char *argv[])
{
	int ret = EXIT_FAILURE;
	try {
		ngila_app app(argc, argv);
		ret = app.run();
	} catch(exception &e) {
		cerror() << e.what() << endl;
	}
	return ret;
}

namespace boost { namespace program_options {
template<>
typed_value<bool>* value() { return bool_switch(); }

template<>
typed_value<bool>* value(bool* v) { return bool_switch(v); }

}}

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
}

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
			cerror() << e.what() << endl;
			throw std::runtime_error("unable to process command line");
	}
}

int ngila_app::run()
{
	if(arg.help || arg.version)
	{
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}

	seq_db mydb;
	for(vector<string>::const_iterator cit = arg_input.begin(); cit != arg_input.end(); ++cit)
	{
		if(!mydb.parse_file(cit->c_str(), true))
		{
			cerror() << "parsing of \'" << cit->c_str() << "\' failed." << endl;
			return EXIT_FAILURE;
		}
	}
	
	return EXIT_SUCCESS;
}

//int main(int argc, char *argv[])
//{
//
//	
//	SeqVec seqs;
//	StringVec names;
//	if(!parse_file(argv[0], names, seqs))
//	{
//		fprintf(stderr, "Parsing \"%s\" failed.\n", argv[0]);
//		return 1;
//	}
//	for_each(seqs.begin(), seqs.end(), seqproc);
//	
//	Sequence seqA, seqB;
//	char msgbuf[255];
//	for(size_t i = 0; i < seqs.size(); ++i)
//	{
//		for(size_t j=i+1; j < seqs.size(); ++j)
//		{
//			double d = align_pair(seqs[i], seqs[j], seqA, seqB);
//			sprintf(msgbuf, "Score = %f", d);
//			print_aln(names[i], seqA, names[j], seqB, msgbuf);
//		}
//	}
//
//
//	return 0;
//}
//
//char letproc(char ch)
//{
//	return (g_bNoCase ? (char)toupper(ch) : ch);
//}
//
//char letproc(const char *cs)
//{
//	return letproc(cs[0]);
//}
//
//double numproc(const char *cs)
//{
//	double d = atof(cs);
//	return (g_bNegate ? -d : d);
//}
//
//void seqproc(string& ss)
//{
//	// Remove Initial Gap
//	string::size_type sy;
//	if(ss[0] == '-')
//	{
//		for(sy = 1; sy < ss.size() && ss[sy] == '-' ; ++sy) { }
//		ss.erase(0, sy);
//	}
//	// Main Body
//	for(string::size_type sz = 0; sz < ss.size(); ++sz)
//	{
//		ss[sz] = letproc(ss[sz]);
//		// Remove Gaps
//		if(ss[sz+1] == '-')
//		{
//			for( sy = sz+2; sy < ss.size() && ss[sy] == '-' ; ++sy)
//				{ }
//			ss.erase(sz+1, sy-sz-1);
//		}
//	}
//}
//
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

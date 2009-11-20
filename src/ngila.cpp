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
#include <sstream>
#include <iomanip>
#include <iostream>

#include <boost/preprocessor.hpp>
#include <boost/foreach.hpp>
#include <boost/config.hpp>

#include "ngila_app.h"
#include "seqdb.h"
#include "matparser.h"
#include "models.h"
#include "align.h"

#define foreach BOOST_FOREACH

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

template<size_t _N>
size_t key_switch(const std::string &ss, const std::string (&key)[_N]) {
	for(size_t i=0;i<_N;++i) {
		if(key[i].find(ss) == 0)
			return i;
	}
	return (size_t)-1;
}

int ngila_app::run()
{
	if(arg.version)	{
		cerr << endl << VERSION_MSG << endl << endl;
		return EXIT_SUCCESS;
	}
	if(arg.help) {
		cerr << endl << VERSION_MSG << endl << endl;
		cerr.precision(5);
		cerr << desc << endl;
		return EXIT_SUCCESS;
	}
	if(arg.quiet)
		cerr.clear(ios::failbit);
	
	seq_db mydb(arg.remove_gaps);
	for(vector<string>::const_iterator cit = arg.input.begin(); cit != arg.input.end(); ++cit) {
		if(!mydb.parse_file(cit->c_str(), true, arg.case_insensitivity)) {
			CERROR("parsing of \'" << cit->c_str() << "\' failed.");
			return EXIT_FAILURE;
		}
	}
	
	int direction = 0;
	if(arg.const_align != 0) {
		direction = mydb.unique_sort(seq_db::DIR_ORI | (arg.const_align & 7));
		mydb.transform(direction);
		if(!(arg.const_align & 8)) { // return to original order
			mydb.rearrange(mydb.db().get<id>().begin());
		}
	}

	cost_model *pmod = NULL;
	string model_keys[] = { string("zeta"), string("geo"), string("cost") };
	switch(key_switch(arg.model, model_keys))
	{
	case 0:
		pmod = new zeta_model;
		break;
	case 1:
		pmod = new geo_model;
		break;
	case 2:
		pmod = new cost_model;
		break;
	default:
		CERROR("unknown model \'" << arg.model << "\'.");
		return EXIT_FAILURE;
	};
	
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
	pair_vec pvec;
	string pairs_keys[] = { string("first"), string("all"), string("each") };
	switch(key_switch(arg.pairs, pairs_keys))
	{
	case 0:
		pvec.push_back(make_pair(0,1));
		break;
	case 1:
		for(size_t i=0;i<mydb.size();++i)
			for(size_t j=i+1;j<mydb.size();++j)
				pvec.push_back(make_pair(i,j));
		break;
	case 2:
		for(size_t i=0;i<mydb.size()-1;i+=2)
			pvec.push_back(make_pair(i,i+1));
		break;
	default:
		CERROR("unknown pairs option \'" << arg.pairs << "\'.");
		return EXIT_FAILURE;
	};
	
	string format_keys[] = { string("aln"), string("fasta") };
	int out_format = 0;
	if(!arg.output.empty()) {
		string ssFormat;
		string::size_type pos = arg.output.find_first_of(':');
#ifdef BOOST_WINDOWS
		if(pos != string::npos && pos != 1) {
#else
		if(pos != string::npos) { // format:file
#endif		
			ssFormat = arg.output.substr(0, pos);
			arg.output.erase(0, pos+1);
		} else { // file.format
			pos = arg.output.find_last_of('.');
			if(pos != string::npos) 
				ssFormat = arg.output.substr(pos+1);
		}
		if(!ssFormat.empty()) {
			out_format = key_switch(ssFormat, format_keys);
			if(out_format == -1) {
				CERROR("unknown output format \'" << ssFormat << "\'.");
				return EXIT_FAILURE;
			}
		}
	}
	
	ofstream fout;
	if(!(arg.output.empty() || arg.output == "-")) {
		fout.open(arg.output.c_str(), ios_base::out|ios_base::trunc);
		if(!fout.is_open()) {
			CERROR("unable to open output file \'" << arg.output << "\'.");
			return EXIT_FAILURE;
		}
	}
	
	for(pair_vec::const_iterator cit = pvec.begin(); cit != pvec.end(); ++cit)
	{
		if(cit != pvec.begin())
			cout << "//" << endl;
		// swap a and b so that a's hashed position is lower
		size_t a = cit->first;
		size_t b = cit->second;
		bool swapped = false;
		if(!(arg.const_align & 8) && mydb.db().project<hashid>(mydb.db().begin()+a) > 
			mydb.db().project<hashid>(mydb.db().begin()+b) ) {
			swap(a,b);
			swapped = true;
		}

		alignment aln(mydb[a], mydb[b]);
		double dcost = alner.align(aln);
			dcost += pmod->offset(mydb[a].dna,
		                      mydb[b].dna);
		
		if(fout.is_open()) {
			aln.print(fout, out_format, dcost,
				((arg.const_align & 16) ? 0 : direction), swapped);
		} else {
			aln.print(cout, out_format, dcost,
				((arg.const_align & 16) ? 0 : direction), swapped);
		}
	}
	
	return EXIT_SUCCESS;
}

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

#include "ngila.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include <boost/preprocessor.hpp>
#include <boost/foreach.hpp>
#include <boost/config.hpp>
//BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <boost/algorithm/string/predicate.hpp>

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
typed_value<bool>* value(bool* v) {
	return bool_switch(v);
}
}}

ngila_app::ngila_app(int argc, char* argv[]) : desc("Allowed Options")
{
	try {
		desc.add_options()
			#define XCMD(lname, sname, desc, type, def) ( \
				_S(lname) _IFD(sname, "," BOOST_PP_STRINGIZE sname), \
				po::value< type >(&arg._V(lname))->default_value(def), \
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

template<std::size_t _N>
std::size_t key_switch(const std::string &ss, const std::string (&key)[_N]) {
	using boost::algorithm::starts_with;
	for(std::size_t i=0;i<_N;++i) {
		if(starts_with(key[i], ss))
			return i;
	}
	return (std::size_t)-1;
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
	seq_db::size_type table_size = mydb.size();
	
	switch(key_switch(arg.pairs, pairs_keys)) {
	case 0:
		pvec.push_back(make_pair(0,1));
		table_size = 2;
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
	
	string format_keys[] = {
		string("aln"), string("fasta"),
		string("dist"), string("dist-c"),
		string("dist-d"), string("dist-i")
	};
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
	ostream &myout = fout.is_open() ? static_cast<ostream&>(fout)
	                                : static_cast<ostream&>(cout);
	typedef boost::multi_array<float, 2> dist_mat;
	dist_mat dist_table;
	const bool do_dist = (out_format >= 2);
	if(do_dist) {
		dist_table.resize(boost::extents[table_size][table_size]);
		fill(dist_table.data(), dist_table.data()+dist_table.num_elements(), numeric_limits<float>::quiet_NaN());
		float diag = (out_format == 5) ? 1.0f : 0.0f;
		for(seq_db::size_type i = 0; i < table_size; ++i)
			dist_table[i][i] = diag;
	}
	
	for(pair_vec::const_iterator cit = pvec.begin(); cit != pvec.end(); ++cit) {
		if(!do_dist && cit != pvec.begin())
			myout << "//" << endl;
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
		dcost += pmod->offset(mydb[a].dna, mydb[b].dna);
		double dident = aln.identity();
		
		if(!do_dist) {
			aln.print(myout, out_format, dcost,
				((arg.const_align & 16) ? 0 : direction), swapped);
		} else {
			double upper,lower;
			if(a > b)
				swap(a,b);
			switch(out_format) {
			default:
			case 2: //dist
				lower = 1.0-dident;
				upper = dcost;
				break;
			case 3: //dist-c
				upper = lower = dcost;
				break;
			case 4: //dist-d
				upper = lower = 1.0-dident;
				break;
			case 5: //dist-i
				upper = lower = dident;
				break;
			}
			dist_table[a][b] = (float)upper;
			dist_table[b][a] = (float)lower;
		}
	}
	if(do_dist) {
		myout << setprecision(10);
		if(arg.header) {
			for(seq_db::size_type i = 0; i < table_size-1; ++i)
				myout << mydb[i].name << "\t";
			myout << mydb[table_size-1].name << endl;
		}
		for(seq_db::size_type i = 0; i < table_size; ++i) {
			for(seq_db::size_type j = 0; j < table_size-1; ++j)
				myout << dist_table[i][j] << "\t";
			myout << dist_table[i][table_size-1] << endl;
		}
	}
	
	if(arg.desktop) {
		cerr << "\nPress any key to continue." << endl;
		char x;
		cin.get(x);
	}
	
	return EXIT_SUCCESS;
}

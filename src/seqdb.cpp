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
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include "seqdb.h"
#include "seqparser.h"

#define foreach BOOST_FOREACH

using namespace std;
using namespace boost;
using namespace boost::multi_index;


bool seq_db::parse_file(const char *csfile, bool bappend, bool bi)
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
		cont.modify(it, bind(&seq_data::sanitize, _1, bi));
	
	return true;
}

struct hash_data {
	std::size_t hashed;
	int dir;
	typedef reference_wrapper<const seq_data> wrap;
	wrap ref;
	
	hash_data(std::size_t h, int d, const seq_data &s) :
		hashed(h), dir(d), ref(s) {
	}
	operator const seq_data&() const {
		return ref;
	}
};

struct full_key : composite_key< hash_data,
	member< hash_data, int, &hash_data::dir>,
	member< hash_data, size_t, &hash_data::hashed>
>{};

typedef multi_index_container< hash_data, indexed_by<
	ordered_non_unique<
		member< hash_data, size_t, &hash_data::hashed>, greater<size_t> >,
	ordered_non_unique<
		full_key, composite_key_result_greater<full_key::result_type> >
	> > hash_container;
	
typedef hash_container::nth_index<1>::type hc_by_both;

int seq_db::unique_sort(int directions) {
	seq_data seq("a", "aa");
	hash_data hd(0, 0, seq);
	const seq_data &ww = (hd);

	hash_container hc;
	foreach(const seq_data & sd, cont) {
		if(directions & DIR_ORI) // hash sequence
			hc.insert(hash_data(hash_range(sd.dna.begin(), sd.dna.end()), DIR_ORI, sd));
		if(directions & DIR_REV) // hash reverse sequence
			hc.insert(hash_data(hash_range(sd.dna.rbegin(), sd.dna.rend()), DIR_REV, sd));
		if(directions & DIR_COM) // hash complement sequence
			hc.insert(hash_data(hash_range(complementer(sd.dna.begin()),
				complementer(sd.dna.end())), DIR_COM, sd));
		if(directions & DIR_RVC) // hash reverse complement
			hc.insert(hash_data(hash_range(complementer(sd.dna.rbegin()),
				complementer(sd.dna.rend())), DIR_RVC, sd));
	}
	size_t hh = hc.begin()->hashed;
	int dd = hc.begin()->dir, mask = -1; 
	foreach(const hash_data &h, hc) {
		if(h.hashed == hh) {
			dd |= h.dir;
			continue;
		}
		dd &= mask;
		switch(dd) {
		case DIR_ORI: case DIR_REV: case DIR_COM: case DIR_RVC:
			break;
		default:
			mask = dd;
		case 0:
			hh = h.hashed;
			dd = h.dir;
			continue;
		}
		break;
	}
	// dd now represents the best direction(s)
	// use an array to break any ties 8 > 1 > 2 > 4
	static const int ties[] = { 8,1,2,1,4,1,2,1,8,8,8,8,8,8,8,8 };
	dd = ties[dd&15];
	// sort
	cont.rearrange(hc.get<1>().find(dd));
	// store this order
	cont.get<hashid>().rearrange(cont.begin());
	return dd;
}

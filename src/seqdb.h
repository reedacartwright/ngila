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

#ifndef SEQDB_H
#define SEQDB_H

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

struct seq_data {
	std::string name;
	std::string dna;
	
	seq_data(const std::string &n, const std::string &s) : name(n), dna(s) {
	}
	
	// append to dna
	seq_data& operator+=(const seq_data& a) {
		dna += a.dna;
		return *this;
	}
	
	// remove characters in dna that belong to g
	seq_data& operator-=(const std::string& g) {
		std::string::size_type (std::string::*gf)(std::string::value_type, std::string::size_type) const
			 = &std::string::find;
		dna.erase( remove_if(dna.begin(), dna.end(),
			boost::bind(gf, g, _1, 0) != std::string::npos),
			dna.end());
	}
	
	void sanitize(bool bi=false) {
		std::replace_if(dna.begin(), dna.end(), std::ptr_fun(::isspace), '_');
		if(bi)
			std::transform(name.begin(), name.end(), name.begin(), std::ptr_fun(::toupper));		
	}
};

//tags
struct name {};

class seq_db {
public:
	typedef boost::multi_index_container< seq_data, boost::multi_index::indexed_by<
		boost::multi_index::random_access<>,
		boost::multi_index::ordered_unique<
			boost::multi_index::tag<name>,
			boost::multi_index::member< seq_data, std::string, &seq_data::name> >
		> > container;
	typedef container::size_type size_type;
	
	inline void add(const seq_data &s)
	{
		// check to see if "name" already exists
		// add to name if it does, otherwise push new sequence
		std::pair<container::iterator,bool> res = cont.push_back(s);
		using boost::lambda::_1;
		if(!res.second)
			cont.modify(res.first, _1 += s);
		// remove gaps
		cont.modify(res.first, _1 -= ss_gaps);
	}
	
	//inline void append(size_type idx, const sequence& s)
	//{
	//	data_vec[idx].second.append(s);
	//}
	
	inline size_type size() const { return cont.size(); }
	
	inline const seq_data& operator[](size_type sz) const {
		return cont[sz];
	}
	
	//inline const data_vec_type& data() const { return data_vec; }
			
	inline void clear() {
		cont.clear();
	}
	
	bool parse_file(const char *csfile, bool bappend=false, bool bi=false);
	
	seq_db() : ss_gaps("-") { }
	seq_db(const std::string &g) : ss_gaps(g) { }
	
protected:
	std::string ss_gaps;
	
	container cont;
};
#endif

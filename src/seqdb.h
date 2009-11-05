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
#include <boost/functional/hash.hpp>

template<class _Iterator>
class _complementer : public _Iterator {
public:
	typedef _Iterator base_type;
	typedef typename base_type::value_type value_type;

	value_type operator*() const {
		static char dnac[] = "\0\x01\x02\x03\x04\x05\x06\a\b\t\n\x0b\f\r\x0e\x0f"
			"\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f"
			" !\"#$%&'()*+,-./0123456789:;<=>?@"
			"TVGHEFCDIJMLKNOPQYSAABWXRZ" "[\\]^_`"
			"tvghefcdijmlknopqysaabwxrz" "{|}~\x7f";
		return value_type(dnac[**static_cast<const base_type*>(this)]);
	}
	
	_complementer(base_type b) : base_type(b) { }
};

template<class _Iterator>
_complementer<_Iterator> complementer(_Iterator it) {
	return _complementer<_Iterator>(it);
}

struct seq_data {
	std::string name;
	std::string dna;
	std::size_t hashed;
	char dir;
	
	enum {
		DIR_ORI = 0, DIR_REV = 1, DIR_COM = 2, DIR_RVC = 4, DIR_ALL = 7, _DIR_ORI = 8
	};
	
	seq_data(const std::string &n, const std::string &s) : name(n), dna(s), hashed(0) {
	}
	
	// append to dna
	inline seq_data& operator+=(const seq_data& a) {
		dna += a.dna;
		return *this;
	}
	
	// remove characters in dna that belong to g
	inline seq_data& operator-=(const std::string& g) {
		std::string::size_type (std::string::*gf)(std::string::value_type, std::string::size_type) const
			 = &std::string::find;
		dna.erase( remove_if(dna.begin(), dna.end(),
			boost::bind(gf, g, _1, 0) != std::string::npos),
			dna.end());
	}
	
	inline void finalize(bool bi=false, int directions=0) {
		std::replace_if(name.begin(), name.end(), std::ptr_fun(::isspace), '_');
		if(bi)
			std::transform(dna.begin(), dna.end(), dna.begin(), std::ptr_fun(::toupper));
		// hash sequence
		dir = _DIR_ORI;
		hashed = boost::hash_range(dna.begin(), dna.end());
		// hash reverse sequence
		if(directions & DIR_REV) {
			std::size_t h = boost::hash_range(dna.rbegin(), dna.rend());
			if(h > hashed) {
				hashed = h;
				dir = DIR_REV;
			} else if(h == hashed) {
				dir |= DIR_REV;
			}
		}
		// hash complement sequence
		if(directions & DIR_COM) {
			std::size_t h = boost::hash_range(
				complementer(dna.begin()), complementer(dna.end()));
			if(h > hashed) {
				hashed = h;
				dir = DIR_COM;
			} else if(h == hashed) {
				dir |= DIR_COM;
			}
		}
		// hash reverse complement
		if(directions & DIR_RVC) {
			std::size_t h = boost::hash_range(
				complementer(dna.rbegin()), complementer(dna.rend()));
			if(h > hashed) {
				hashed = h;
				dir = DIR_RVC;
			} else if(h == hashed) {
				dir |= DIR_RVC;
			}
		}
	}
};

//tags
struct name {};
struct id {};
struct hashid {};

class seq_db {
public:
	typedef boost::multi_index_container< seq_data, boost::multi_index::indexed_by<
		boost::multi_index::random_access<>,
		boost::multi_index::random_access<boost::multi_index::tag<id> >,
		boost::multi_index::ordered_unique<
			boost::multi_index::tag<name>,
			boost::multi_index::member< seq_data, std::string, &seq_data::name> >,
		boost::multi_index::ordered_non_unique<
			boost::multi_index::tag<hashid>,
			boost::multi_index::member< seq_data, std::size_t, &seq_data::hashed>,
			std::greater<std::size_t> >
	> > container;
	typedef container::size_type size_type;
	
	inline void add(const seq_data &s) {
		// check to see if "name" already exists
		// add to name if it does, otherwise push new sequence
		std::pair<container::iterator,bool> res = cont.push_back(s);
		using boost::lambda::_1;
		if(!res.second)
			cont.modify(res.first, _1 += s);
		// remove gaps
		cont.modify(res.first, _1 -= ss_gaps);
	}
		
	inline size_type size() const { return cont.size(); }
	
	inline const seq_data& operator[](size_type sz) const {
		return cont[sz];
	}
				
	inline void clear() {
		cont.clear();
	}
	
	bool parse_file(const char *csfile, bool bappend=false, bool bi=false, int dirs=0);
	
	seq_db() : ss_gaps("-") { }
	seq_db(const std::string &g) : ss_gaps(g) { }
	
protected:
	std::string ss_gaps;
	
	container cont;
};


#endif

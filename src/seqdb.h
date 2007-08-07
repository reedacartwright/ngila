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
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>

class seq_db
{
public:
	typedef std::string sequence;
	typedef std::string name;
	typedef std::pair<name, sequence> value_type; 
	typedef std::vector<value_type>::size_type size_type;
	typedef std::map<name, size_type> name_map;
	
	inline bool add(const name& id, const sequence& s)
	{
		sequence t;
		std::remove_copy_if(s.begin(), s.end(), std::back_inserter(t),
			boost::bind(boost::mem_fn<sequence::size_type,sequence, 
			sequence::value_type, sequence::size_type>(&sequence::find),
			ss_gaps, _1, 0) != sequence::npos
		);
// 		std::remove_copy_if(s.begin(), s.end(), std::back_inserter(t),
// 		                    boost::bind(strchr, ss_gaps.c_str(), _1) != (char*)NULL
// 		);
// 		std::remove_copy_if(s.begin(), s.end(), std::back_inserter(t),
// 		                    boost::lambda::_1 == '-');

		
		std::pair<name_map::iterator, bool> p = data_map.insert(make_pair(id, size()));
		
		if(p.second)
			data_vec.push_back(make_pair(id, t));
		else
			data_vec[p.first->second].second.append(t);
		return p.second;
	}
	inline void append(size_type idx, const sequence& s)
	{
		data_vec[idx].second.append(s);
	}
	
	inline size_type size() const { return data_vec.size(); }
	
	inline const value_type& operator[](size_type sz) const
	{
		return data_vec[sz];
	}
			
	inline void clear()
	{
		data_vec.clear();
		data_map.clear();
	}
	
	bool parse_file(const char *csfile, bool bappend=false);
	
	seq_db() : ss_gaps("-") { }
	seq_db(const std::string &g) : ss_gaps(g) { }
	
protected:
	std::vector<value_type> data_vec;
	name_map data_map;
	std::string ss_gaps;
};

#endif

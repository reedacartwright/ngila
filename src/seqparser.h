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

#include <iostream>
#include <stack>

#include "seqdb.h"

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_core.hpp>
using namespace boost::spirit::classic;

struct push_string
{
	push_string(std::stack<std::string>& st) : my_stack(st) { }
	
	template<typename It> void operator() ( It first, It last) const
	{
		my_stack.push(std::string(first, last));
	}
		
	std::stack<std::string>& my_stack;
};

inline void trim_string(std::string &ss, const char *ch = " \n\t\r\v\f")
{
	std::string::size_type sz = ss.find_last_not_of(ch);
	if(sz != std::string::npos) ss.erase(sz+1); 
	sz = ss.find_first_not_of(ch);
	if(sz != std::string::npos) ss.erase(0, sz);
}

inline void sanitize_string(std::string &ss, const char *ch = " \n\t\r\v\f")
{
	std::string::size_type sz = ss.find_first_of(ch);
	while(sz != std::string::npos)
	{
		ss.erase(sz, 1);
		sz = ss.find_first_of(ch, sz);
	}
}

struct pop_sequence
{
	pop_sequence(std::stack<std::string>& st, seq_db& db) : my_stack(st), rdb(db) { }
	
	template<typename It> void operator() ( It first, It last) const
	{
		std::string seq(my_stack.top());
		my_stack.pop();
		std::string name(my_stack.top());
		my_stack.pop();
		trim_string(name);
		sanitize_string(seq);
		
		rdb.add(seq_data(name, seq));
	}
		
	std::stack<std::string>& my_stack;
	seq_db &rdb;
};

struct seq_grammar : public grammar<seq_grammar>
{
	seq_grammar(std::stack<std::string>& st, seq_db& db ) : string_stack(st), rdb(db) {}
	
	template <typename ScannerT> struct definition
	{
		definition(seq_grammar const& self)
		{
			file_format =
					(*space_p) >>
					fasta_format
//				|	clustal
//				|	phylip
				;
			fasta_format = 
					*fasta_seq
				;
			fasta_seq =
				(	fasta_seq_head
				>>	fasta_seq_body
				)	[pop_sequence(self.string_stack, self.rdb)]
					;
			fasta_seq_head = 
					ch_p('>')
				>>	(+(graph_p|blank_p))[push_string(self.string_stack)]
				>>	eol_p
				;
			fasta_seq_body =
					(+(~ch_p('>')))[push_string(self.string_stack)]
				;
			
		}

		rule<ScannerT> file_format;
		rule<ScannerT> fasta_format;
		rule<ScannerT> fasta_seq;
		rule<ScannerT> fasta_seq_head;
		rule<ScannerT> fasta_seq_body;
		rule<ScannerT> const& start() const { return file_format; }
	};
	
	std::stack<std::string> &string_stack;
	seq_db &rdb;
};

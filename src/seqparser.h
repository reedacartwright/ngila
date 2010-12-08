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

struct add_sequence
{
	add_sequence(seq_db::size_type &pos, seq_db& db) : rpos(pos), rdb(db) { }
	
	template<typename It> void operator() ( It first, It last) const
	{
		std::string seq(first,last);
		sanitize_string(seq);
		seq_db::size_type pos = rpos++ % rdb.size();
		rdb.add(seq_data(rdb[pos].name, seq));
	}
	
	seq_db::size_type &rpos;
	seq_db &rdb;
};


struct seq_grammar : public grammar<seq_grammar> {
	seq_grammar(std::stack<std::string>& st, seq_db& db, seq_db::size_type &p) : string_stack(st), rdb(db), pos(p) {}
	
	template <typename ScannerT> struct definition {
		definition(seq_grammar const& self) {
			self.pos = 0;
			
			file_format = (*space_p) >> (phylip_format|aln_format|fasta_format)
				>> (*space_p);
			blank_line = *blank_p >> eol_p;

			fasta_format = *fasta_seq;
			fasta_seq =	(fasta_seq_head >> fasta_seq_body)
				[pop_sequence(self.string_stack, self.rdb)];
			fasta_seq_head = ch_p('>')
				>>(+(graph_p|blank_p))[push_string(self.string_stack)]
				>>	eol_p
				;
			fasta_seq_body = (+(~ch_p('>')))[push_string(self.string_stack)];

			aln_format = aln_head >> +aln_line;
			aln_head = str_p("CLUSTAL") >> *(graph_p|blank_p) >> eol_p;
			aln_line = (aln_seq | aln_special | *blank_p) >> eol_p;
			aln_seq = (aln_seq_name >> +blank_p >> aln_seq_body >> *blank_p
				>> *digit_p)
				[pop_sequence(self.string_stack, self.rdb)];
			aln_seq_name = (+graph_p)[push_string(self.string_stack)];
			aln_seq_body = (+graph_p)[push_string(self.string_stack)];
			aln_special = +blank_p >> punct_p >> *(blank_p|punct_p);

			phylip_format =	phylip_head >> phylip_block1
				>> *(blank_line >> phylip_block);
			phylip_head = uint_p >> +(blank_p) >> uint_p >> (*blank_p >> eol_p);
			phylip_block1 = *(phylip_seq1 >> eol_p);
			phylip_seq1 = (phylip_seq_name >> phylip_seq_seq)
				[pop_sequence(self.string_stack, self.rdb)];
			phylip_seq_name = (repeat_p(10)[graph_p|blank_p])[push_string(self.string_stack)];
			phylip_seq_seq = (+(graph_p|blank_p))[push_string(self.string_stack)];
			phylip_block = *(*blank_p >> phylip_seq >> eol_p);
			phylip_seq = (+(graph_p|blank_p))[add_sequence(self.pos, self.rdb)];
			
		}

		rule<ScannerT> file_format, blank_line;
		rule<ScannerT> fasta_format, fasta_seq, fasta_seq_head, fasta_seq_body;
		rule<ScannerT> aln_format, aln_head, aln_line, aln_seq, aln_seq_name,
			aln_seq_body, aln_special;
		rule<ScannerT> phylip_format, phylip_head, phylip_block1, phylip_block,
			phylip_seq1, phylip_seq, phylip_seq_name, phylip_seq_seq;
		
		rule<ScannerT> const& start() const { return file_format; }
	};
	
	std::stack<std::string> &string_stack;
	seq_db::size_type &pos;
	seq_db &rdb;
};

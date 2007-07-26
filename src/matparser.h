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

#include <boost/spirit/core.hpp>
#include <boost/spirit.hpp>

using namespace boost::spirit;

struct mat_grammar : public grammar<mat_grammar>
{
	mat_grammar() : {}
	
	template <typename ScannerT> struct definition
	{
		definition(mat_grammar const& self)
		{
			matrix = header >> body;
			header = *graph_p;
			body = *line;
			line = word_p >> *real_p;

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
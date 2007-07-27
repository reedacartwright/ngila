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
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/swap_actor.hpp>
#include <boost/spirit/actor/clear_actor.hpp>
#include <boost/spirit.hpp>

typedef double sub_matrix[128][128];

bool parse_matrix(const char *cs, sub_matrix &rsm);

struct mat_work
{
	typedef std::vector<double> row_type;
	typedef std::vector<row_type> mat_type;
	typedef std::string label_type;

	label_type labels;
	mat_type data;
	row_type rtemp;
	
	void clear()
	{
		labels.clear();
		data.clear();
		rtemp.clear();
	}

	bool process(sub_matrix &m);
};

using namespace boost::spirit;

struct mat_grammar : public grammar<mat_grammar>
{
	mat_grammar(mat_work& w) : work(w) { work.clear(); }
	
	template <typename ScannerT> struct definition
	{
		definition(mat_grammar const& self)
		{
			matrix = *(anychar_p) >> end_p;
			header = (*graph_p)[assign_a(self.work.labels)] >> eol_p;
			body = *line[push_back_a(self.work.data, self.work.rtemp)]
				[clear_a(self.work.rtemp)];
			line = !graph_p >> *real_p[push_back_a(self.work.rtemp)] >> eol_p;
		}

		rule<ScannerT> matrix;
		rule<ScannerT> header;
		rule<ScannerT> body;
		rule<ScannerT> line;
		rule<ScannerT> const& start() const { return matrix; }
	};
	
	mat_work &work;
};
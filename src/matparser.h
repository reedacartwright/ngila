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

#include <boost/spirit.hpp>
#include <boost/spirit/core.hpp>
#include <boost/spirit/core/parser.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/clear_actor.hpp>

typedef double sub_matrix[128][128];

bool parse_matrix(const char *cs, sub_matrix &rsm);

struct mat_work
{
	typedef std::vector<char> labels_type;
	typedef std::vector<double> row_type;
	typedef std::vector<row_type> data_type;
	typedef sub_matrix matrix;

	labels_type labels;
	data_type data;
	row_type rt;

	inline void clear()
	{
		labels.clear();
		data.clear();
		rt.clear();
	}

	bool process( matrix &m) const;
};

using namespace boost::spirit;

struct mat_grammar : public grammar<mat_grammar>
{
	mat_grammar(mat_work &w) : work(w) { work.clear(); }
	
	mat_work &work;

	template <typename ScannerT> struct definition
	{
		definition(mat_grammar const& self)
		{
			matrix = *eol_p >> header >> body >> end_p;
			header = +(graph_p[push_back_a(self.work.labels)]);
			body = +line >> !eol_p;
			line = eol_p >> !(graph_p - real_p) >>
				(+(real_p[push_back_a(self.work.rt)]))
				[push_back_a(self.work.data, self.work.rt)]
				[clear_a(self.work.rt)]
				;
		}

		rule<ScannerT> matrix;
		rule<ScannerT> header;
		rule<ScannerT> body;
		rule<ScannerT> line;
		rule<ScannerT> const& start() const { return matrix; }
	};
};


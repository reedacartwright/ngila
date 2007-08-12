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

#include <boost/preprocessor.hpp>

#ifndef XMACROS_HELPERS
#define XMACROS_HELPERS

/****************************************************************************
 *    X-Helpers List                                                        *
 ****************************************************************************/

// The _JS macro cats a seq 'seq' with separator 'sep'

#define _JS_OP(s, data, elem) BOOST_PP_CAT(data, elem)

#define _JS(sep, seq) BOOST_PP_IF( BOOST_PP_EQUAL(BOOST_PP_SEQ_SIZE(seq),1), \
	_JS_1, _JS_2)(sep,seq)

#define _JS_1(sep, seq) BOOST_PP_SEQ_HEAD(seq)

#define _JS_2(sep, seq) BOOST_PP_SEQ_CAT(( BOOST_PP_SEQ_HEAD(seq) ) \
	BOOST_PP_SEQ_TRANSFORM(_JS_OP, sep, BOOST_PP_SEQ_TAIL(seq)) \
)

// The _SS macro is similiar to _JS except that it stringizes everything

#define _SS_OP(r, data, elem) data BOOST_PP_STRINGIZE(elem)

#define _SS(sep, seq) BOOST_PP_STRINGIZE(BOOST_PP_SEQ_HEAD(seq)) \
	BOOST_PP_SEQ_FOR_EACH(_SS_OP, sep, BOOST_PP_SEQ_TAIL(seq))

// Tests whether a seq is defined.

#define _SD(seq) BOOST_PP_GREATER(BOOST_PP_SEQ_SIZE(seq),1)

#else

/***************************************************************************
 *    Cleanup                                                              *
 ***************************************************************************/

#undef XMACROS_HELPERS
#undef _JS_OP
#undef _JS
#undef _SS
#undef _SD

#endif


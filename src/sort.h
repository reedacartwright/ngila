#pragma once
#ifndef NGILA_SORT_H
#define NGILA_SORT_H

#include <algorithm>
#include <vector>
#include <string>

template<class Base>
struct seq_sort_element {
	seq_sort_element(const Base& s) : seq(s), direction(0), order(-1) { }
	Base seq;
	std::size_t order;
	std::size_t orig_order;
	std::size_t direction;
};

template<class _Seq>
class seq_sort : public std::vector< seq_sort_element<_Seq> >  {
public:
	typedef seq_sort_element<_Seq> element_type;
	typedef std::vector< seq_sort_element<_Seq> > base_type;

	enum {
		REV = 1,
		COM = 2,
		RVC = 4,
		ALL = 7
	};
	
	seq_sort() { }
	
	template<It>
	seq_sort(It first, It last) : base_type(first, last) {
		base_type::size_type o = 0;
		for(base_type::iterator it = begin(); it != end(); ++it) {
			it->orig_order = o++;
		}
	}

	int seq_sort(int dirs) {
		
	}
};

inline char dna_complement(char ch) {
	static char dnac[] = "\0\x01\x02\x03\x04\x05\x06\a\b\t\n\x0b\f\r\x0e\x0f"
		"\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f"
		" !\"#\$%&'()*+,-./0123456789:;<=>?@"
		"TVGHEFCDIJMLKNOPQYSAABWXRZ" "[\\]^_`"
		"tvghefcdijmlknopqysaabwxrz" "{|}~\x7f";
	return dnac[ch];
}

template<It>
void seq_transform(It first, It last, int dir) {
	if(dir == 0) {
		/*noop*/;
	} else if(dir == 1) {
		std::reverse(first,last);
	} else if(dir == 2) {
		std::transform(first, last, first, dna_complement);
	} else if(dir == 4 || dir == 3) {
		std::reverse(first, last);
		std::transform(first, last, first, dna_complement);		
	}
}

ngila_sort sorter(a,b, ngila_sort::ALL);
sorter.val[0];
sorter.val[1];


#endif // NGILA_SORT_H


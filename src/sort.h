#pragma once
#ifndef NGILA_SORT_H
#define NGILA_SORT_H

#include <algorithm>
#include <vector>
#include <string>
#include <boost/functional/hash.hpp>

enum {
	DIR_ORI = 0,
	DIR_REV = 1,
	DIR_COM = 2,
	DIR_RVC = 4,
	DIR_ALL = 7
};

inline char dna_complement(char ch) {
	static char dnac[] = "\0\x01\x02\x03\x04\x05\x06\a\b\t\n\x0b\f\r\x0e\x0f"
		"\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f"
		" !\"#$%&'()*+,-./0123456789:;<=>?@"
		"TVGHEFCDIJMLKNOPQYSAABWXRZ" "[\\]^_`"
		"tvghefcdijmlknopqysaabwxrz" "{|}~\x7f";
	return dnac[ch];
}

template<class It>
inline void seq_transform(It first, It last, int dir) {
	switch(dir&7) {
	case DIR_REV:
		std::reverse(first,last);
		break;
	case DIR_RVC:
	case (DIR_COM|DIR_REV):
		std::reverse(first, last);
	case DIR_COM:
		std::transform(first, last, first, dna_complement);
	case DIR_ORI:
	default:
		break;
	}
}

template<class Base>
struct seq_sort_element {
	seq_sort_element(const Base& s) : seq(s), direction(0) { }
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
	typedef _Seq value_type;
	
	seq_sort() { }
	
	template<class It>
	seq_sort(It first, It last) : base_type(first, last) {
		typename base_type::size_type o = 0;
		for(typename base_type::iterator it = this->begin(); it != this->end(); ++it) {
			it->order = it->orig_order = o++;
		}
	}
	
	struct hash_data {
		std::size_t hash_value;
		std::size_t orig_order;
		std::size_t direction;
		hash_data(std::size_t a, std::size_t b, std::size_t c) :
			hash_value(a), orig_order(b), direction(c) { }
			
		bool operator<(const hash_data &h) const {
			return hash_value < h.hash_value;
		}
	};
	
	// Dirs is a flag that specifies which directions to try
	int sort_seqs(int dirs) {
		// create hashing object
		boost::hash<value_type> hasher;
		// store is a buffer of hashes, which will be used to determine 
		// unique ordering and direction
		std::vector<hash_data> store;
		std::size_t o = 0;
		value_type ss;
		
		// fill store with hashed values
		for(typename base_type::iterator it = this->begin(); it != this->end(); ++it,++o) {
			// original
			store.push_back(hash_data(hasher(it->seq), o, DIR_ORI));
			dirs &= 7;
			if(dirs == DIR_ORI)
				continue;
			ss = it->seq;
			if(dirs & DIR_COM == DIR_ORI) { // no complement sequence requested
				seq_transform(ss.begin(), ss.end(), DIR_REV);
				if(dirs & DIR_REV)
					store.push_back(hash_data(hasher(ss), o, DIR_REV));
				if(dirs & DIR_RVC) {
					seq_transform(ss.begin(), ss.end(), DIR_COM);
					store.push_back(hash_data(hasher(ss), o, DIR_RVC));
				}
			} else { // complement sequence is requested
				seq_transform(ss.begin(), ss.end(), DIR_COM);
				store.push_back(hash_data(hasher(ss), o, DIR_COM));
				if(dirs != DIR_COM) { // other sequences requested
					seq_transform(ss.begin(), ss.end(), DIR_REV);
					if(dirs & DIR_RVC)
						store.push_back(hash_data(hasher(ss), o, DIR_RVC));
					if(dirs & DIR_REV) {
						seq_transform(ss.begin(), ss.end(), DIR_COM);
						store.push_back(hash_data(hasher(ss), o, DIR_REV));
					}
				}
					
			}
		}
		// sort the hashed values
		std::sort(store.begin(), store.end());
		
		//lets see the look
		for(typename std::vector<hash_data>::iterator it = store.begin(); it != store.end(); ++it) {
			std::cerr << it->hash_value << " "
				<< it->orig_order << " "
				<< it->direction << std::endl;
		}
	}
};

#endif // NGILA_SORT_H


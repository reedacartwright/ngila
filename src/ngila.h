#ifndef NGILA_H
#define NGILA_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <vector>
#include <algorithm>
#include <map>

// Typedefs
typedef std::string Sequence;
typedef std::vector<Sequence> SeqVec;
typedef std::vector<std::string> StringVec;
typedef std::map<std::string, Sequence> SeqMap;

// External variables
extern double mCost[128][128]; //substitution scoring
extern double dB,dA,dC;

// Static Variables
static const char chGap = '-';

// Functions
double align_pair(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD); // align A and B, placing results in C and D
double align_pair_x(const Sequence& seqA, const Sequence& seqB, Sequence& seqC, Sequence& seqD);
char letproc(const char *cs);
double numproc(const char *cs);
void seqproc(std::string& ss);

bool parse_matrix(const char* csFile);
bool parse_file(const char* csFile, StringVec &vNames, SeqVec &vSeqs);

// Inline Functions
inline void delete_func(void *p) { delete p; }

// Templates Functions
template<typename val> val min3(val a, val b, val c)
{
	return std::min(a,std::min(b,c));
}
template<typename val, typename rep> val min3(val a, val b, val c, rep& r)
{
	if(a <= c)
		if(a <= b)
			{ r = (rep)0; return a; }
		else
			{ r = (rep)1; return b; }
	else if(b <= c)
		{ r = (rep)1; return b; }
	else
		{ r = (rep)2; return c; }
}

// Classes
class SeqDB
{
public:
	typedef SeqVec::size_type Pos;
	typedef std::map<std::string, Pos> NameMap;
	
	const static Pos npos = static_cast<Pos>(-1);

	Pos LookupPos(const std::string& ss) const
	{
		NameMap::const_iterator cit = m_map.find(ss);
		if(cit == m_map.end())
			return npos;
		else
			return cit->second;
	}

	bool Add(const std::string& ss, const Sequence& seq)
	{
		Pos p = LookupPos(ss);
		if(p == npos)
		{
			m_names.push_back(ss);
			m_seqs.push_back(seq);
			m_map[ss] = m_names.size()-1;
		}
		else
		{
			m_seqs[p].append(seq);
		}
		return true;
	}

	bool Add(Pos p, const Sequence& seq)
	{
		if(p >= m_seqs.size())
			return false;
		m_seqs[p].append(seq);
		return true;
	}

	void Clear()
	{
		m_seqs.clear();
		m_names.clear();
		m_map.clear();
	}

	void Transfer(StringVec &vNames, SeqVec &vSeqs)
	{
		std::swap(vNames, m_names);
		std::swap(vSeqs, m_seqs);
		Clear();
	}
	

protected:
	SeqVec m_seqs;
	StringVec m_names;
	NameMap m_map;
};


#endif //NGILA_H

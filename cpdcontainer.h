#ifndef CPDCONTAINER_H_INCLUDED
#define CPDCONTAINER_H_INCLUDED
#include "compound.h"
#include <set>

#include <iostream>
//#define _DEBUG  //already defined in main.cpp

struct AreaCmp{
    bool operator() (const Compound& c1, const Compound& c2) const{
	    return c1.Area > c2.Area;
	}
};

struct MzCmp{
	bool operator() (const Compound& cpd1, const Compound& cpd2) const{
		return cpd1.Mz < cpd2.Mz;
	}
};


typedef multiset<Compound, AreaCmp>  AreaOrderedSet;
typedef multiset<Compound, MzCmp>    MzOrderedSet;
typedef multiset<Compound>::iterator CompoundItr;
typedef tuple<int, vector<Compound> > CpdCluster;


class UnalignedList {
public:
	AreaOrderedSet areaRank;
	MzOrderedSet   mzRank;
public:
	size_t size();
	bool empty();
	void insert(Compound x);
	void remove(const Compound& x);
	void mzSearch(double lb, double ub, CompoundItr& beginItr, CompoundItr& endItr);
	CompoundItr areaSearch(Compound x);
	Compound getTopArea();
};
size_t UnalignedList::size() {
	return areaRank.size();
}

bool UnalignedList::empty() {
	return areaRank.empty();
}

void UnalignedList::insert(Compound x){
	areaRank.insert(x);
	mzRank.insert(x);
}

void UnalignedList::remove(const Compound & x){
	auto areaItr=areaRank.find(x);
	auto mzItr = mzRank.find(x);
	if (areaItr != areaRank.end())
		areaRank.erase(areaItr);
	if (mzItr != mzRank.end())
		mzRank.erase(mzItr);
}

void UnalignedList::mzSearch(double lb, double ub, CompoundItr & beginItr, CompoundItr & endItr){
    Compound lbCpd = Compound(lb, 0, 0);
    Compound ubCpd = Compound(ub, 0, 0);
#ifdef _DEBUG
    std::cout << "    function mzSearch( ) internal\n    {" << endl;
    std::cout << "      Lowerbound Compound:\n";
    std::cout << "        " << lbCpd.toString();
    std::cout << "      Upperbound Compound:\n";
    std::cout << "        " << ubCpd.toString();
    std::cout << "    }" << endl;
#endif

	beginItr = mzRank.lower_bound(lbCpd);
	endItr   = mzRank.upper_bound(ubCpd);
}

CompoundItr UnalignedList::areaSearch(Compound x){
	return areaRank.find(x);
}

Compound UnalignedList::getTopArea(){
	return *(areaRank.begin());
}
#endif // CPDCONTAINER_H_INCLUDED

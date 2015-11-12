#ifndef COMPOUNDTYPE_H_INCLUDED
#define COMPOUNDTYPE_H_INCLUDED
#include<string>
#include<cstdio>

class Compound {
public:
    int dsID;
    std::string dsName;

    int cID;
    std::string cName;

    std::string additionalData;

    double Mz, Rt, Area;
public:
    Compound(double mz,double rt, double area):Mz(mz),Rt(rt),Area(area){}
    Compound(int cid, double mz, double rt, double area):cID(cid), Mz(mz), Rt(rt), Area(area) {}
    Compound(){}
    std::string toString(bool abbr);
};

std::string Compound::toString(bool abbr=false){
    char buf[512] = {0};
    sprintf(buf,"dsID=%d cID=%d m/z=%.10f  Rt=%.10f  Area=%.10f\n",
            dsID, cID, Mz, Rt, Area);
    return std::string(buf);
}

bool CmpByIncreasingRt(Compound cpd1, Compound cpd2) {
	return cpd1.Rt < cpd2.Rt;
}

bool CmpByDsIDArea(Compound cpd1, Compound cpd2){
    if (cpd1.dsID!=cpd2.dsID){
        return cpd1.dsID<cpd2.dsID;
    }else{
        return cpd1.Area>cpd2.Area;
    }
}


#endif // COMPOUNDTYPE_H_INCLUDED

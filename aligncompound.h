#ifndef ALIGNCOMPOUND_H_INCLUDED
#define ALIGNCOMPOUND_H_INCLUDED

#include <vector>
#include <string>

inline string i2s(int n){
    string ret;
    stringstream ss;
    ss<<n;
    ss>>ret;
    return ret;
}


class Peak{
public:
    int dsID,cID;
    Peak(int ds,int c):dsID(ds),cID(c){}
public:
    std::string toString(bool xmlmode);
};
string Peak::toString(bool xmlmode=false){
    /*
    <Peaks>
        <Peak>
            <dsID>0</dsID>
            <cID>1221</cID>
        </Peak>
    </Peaks>
    */
    std::string ret="";
    if (xmlmode==true){
        ret+="<Peak>\n";
        ret+="    <dsID>"+i2s(dsID)+"</dsID>\n";
        ret+="    <cID>" +i2s(cID) +"</cID>\n";
        ret+="</Peak>\n";

    }else{
        ret+="peak: dsID="+i2s(dsID)+" cID="+i2s(cID)+"\n";
    }
    return ret;
}


typedef vector<Peak> PeakContainer;
class AlignedCompound{
public:
    int AlignID;
    PeakContainer Peaks;
public:
    std::string toString(bool xmlmode);
};
std::string AlignedCompound::toString(bool xmlmode=false){
    std::string ret="";
    if (xmlmode==true){
        ret+="<Peaks>\n";
        ret+="<AlignID>"+i2s(AlignID)+"</AlignID>\n";
        for(auto peak : Peaks){
            ret+=peak.toString(true);
        }
        ret+="</Peaks>\n\n";
    }else{
        ret+="AlignID="+i2s(AlignID)+"\n";
        for(auto peak:Peaks){
            ret+="  "+peak.toString();
        }
        ret+="\n\n";
    }
    return ret;
}


typedef vector<AlignedCompound> AlignedList;
#endif // ALIGNCOMPOUND_H_INCLUDED

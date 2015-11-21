#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>

//#define _DEBUG

#include "csvproc.h"
#include "compound.h"
#include "cpdcontainer.h"
#include "aligncompound.h"
using namespace std;

void strvecOutput(const vector<string>& v,
                  bool byCol=false, string sep=" "){
    if (byCol) sep="\n";
    for (auto item:v){
        cout<<item<<sep;
    }
    cout<<endl;
}

inline double s2f(string s){
    double num;
    stringstream ss(s);
    ss>>num;
    return num;
}

#define DATA_INPUT_ERROR 1

int main(int argc, char* argv[]){
	/*
    CsvFile c1("csvtest.csv");
    cout<<c1.table(2,22)<<endl<<"ee"<<endl;
    vector<string> row0=c1.column(22);
    //strvecOutput(row0,true);
    if (c1.table(2,22)=="")
    cout << "Hello world!" << endl;
    system("ls -al");
    */
    double AlignWindowPhase1;
    double AlignWindowPhase2;
    double MassTol;

	vector<CsvFile> dsGroups;
    vector<size_t> datBeginIndices;
    vector<size_t> datEndIndices;

#ifdef _DEBUG
    cout<<"argument cnt:"<<argc<<endl;
    for (int i=0;i<argc;++i){
        cout<<"arg"<<i<<": "<<argv[i]<<endl;
    }
#endif

	/// Get parameters and filename from arguments
	if (argc>=4){
        AlignWindowPhase1=s2f(argv[1]);
        AlignWindowPhase2=s2f(argv[2]);
        MassTol          =s2f(argv[3]);
	}else{
        cerr<<"Data input error: parameter not specified!"<<endl;
        return DATA_INPUT_ERROR;
	}
	for (int i=4;i<argc;++i){
		CsvFile csv(argv[i]);
        dsGroups.push_back(csv);
	}

    /// Define Two Major Containers
    UnalignedList unalignList;
    AlignedList   alignList;

    ///
    for(size_t dsID=0; dsID!=dsGroups.size();++dsID) {
        ///CSV format :
        ///Col:   0  1  2      3             cnt    n        n+1
        ///Title: mz rt Expr1f Expr2f .....  Exprnf isotopes adduct
        ///       #  #  LAB1   LAB1   .....  LAB2   #        #

        /// Get columns containing data, return datBeginIdx& datEndIdx
        CsvFile ds=dsGroups[dsID];
        vector<string> titleRow=ds.row(0);
        size_t datBeginIdx=2, datEndIdx=0;
        for (size_t i=datBeginIdx; i!=titleRow.size();++i){
            if (titleRow[i]=="isotopes"){
                datEndIdx=i;
                break;
            }
        }
        if (datEndIdx==0){
            cerr<<"Data input error: column 'isotopes' not found"<<endl;
            return DATA_INPUT_ERROR;//data format error
        }
        --datEndIdx;
        datBeginIndices.push_back(datBeginIdx);
        datEndIndices.push_back(datEndIdx);

        /// Enumerate every compound
        for (size_t row=2; row<ds.xsize(); ++row){
            Compound cpd;
            cpd.dsID=dsID;
            cpd.cID =(int)row-2;
            cpd.Mz  =s2f(ds.table(row,0));
            cpd.Rt  =s2f(ds.table(row,1));
            double maxArea =s2f(ds.table(row,datBeginIdx));
            for (size_t col=datBeginIdx+1; col<=datEndIdx; ++col){
                maxArea=std::max(maxArea, s2f(ds.table(row,col)));
            }
            cpd.Area=maxArea;
            unalignList.insert(cpd);

        }


    }

    int AlignID = 0;
    int AppearTimes=dsGroups.size();

#ifdef _DEBUG
    cout <<"AppearTime set to "<< AppearTimes << endl;
    cout << "press any key to enter alignment cycling" << endl;
    system("pause");
    int cyccnt = 0;
#endif

	while (!unalignList.empty()) {

#ifdef _DEBUG
        cout << "Iteration times:" << cyccnt++ << "  unalignList.size" << unalignList.size() << endl;
#endif

		Compound cpd = unalignList.getTopArea();
		double AlignRt       = cpd.Rt;
		double AlignMoleMass = cpd.Mz;

#ifdef _DEBUG
        cout << "    primary reference compound:" << endl;
        cout << "      "<< cpd.toString();
#endif

#ifdef _DEBUG
        cout << "    mass seraching range:" << endl;
        cout << "      [" << AlignMoleMass - MassTol << " , " << AlignMoleMass + MassTol << "]" << endl;

#endif

		vector<Compound> alignGroup;
		CompoundItr itrBegin=unalignList.mzRank.end(), itrEnd= unalignList.mzRank.end();
		unalignList.mzSearch(AlignMoleMass - MassTol, AlignMoleMass + MassTol, itrBegin, itrEnd);//find Mass-Valid compounds

#ifdef _DEBUG
        if (itrBegin == unalignList.mzRank.end())
            cout << "    itrBegin Not Found" << endl;
        if (itrEnd == unalignList.mzRank.end())
            cout << "    itrEnd Not Found" << endl;

        cpd = *(itrBegin);
        cout << "    range-left end:" << endl;
        cout << "      " << cpd.toString() << endl;
        cpd = *(itrEnd);
        cout << "    range-right end:" << endl;
        cout << "      " << cpd.toString() << endl;
#endif

		for (auto itr = itrBegin; itr != itrEnd; ++itr) {//find Rt-Valid compounds within Mass-Valid compounds
			auto cpd = *itr;
			if (AlignRt - AlignWindowPhase1*0.5 < cpd.Rt && cpd.Rt < AlignRt + AlignWindowPhase1*0.5)
				alignGroup.push_back(cpd);
		}

		for (auto cpd : alignGroup)
			unalignList.remove(cpd);

		int alignSize = alignGroup.size();
		if ((int)alignSize >= AppearTimes) {//only align the compounds when it is contained in all groups
			sort(alignGroup.begin(), alignGroup.end(), CmpByIncreasingRt);
			AlignRt = alignGroup[alignSize / 2].Rt; //find the compound with median Rt
			AlignMoleMass = alignGroup[alignSize / 2].Mz;

			//enlargement of compounds after the AlignRt and the AlignMoleMass has been revised
			unalignList.mzSearch(AlignMoleMass - MassTol, AlignMoleMass + MassTol, itrBegin, itrEnd);
			for (auto itr = itrBegin; itr != itrEnd; ++itr) {
				auto cpd = *itr;
				if (AlignRt - AlignWindowPhase2*0.5 < cpd.Rt && cpd.Rt < AlignRt + AlignWindowPhase2*0.5)//The Rt bound changes
					alignGroup.push_back(cpd);
			}
		}

		//remove selected compounds
		for (auto cpd : alignGroup)
			unalignList.remove(cpd);

		//wash the aligning result
		if ((int)alignGroup.size() >= AppearTimes) {
			sort(alignGroup.begin(), alignGroup.end(), CmpByDsIDArea);

			PeakContainer pc;

			Compound prevCpd;
			prevCpd.dsID=-1;
            for (auto cpd: alignGroup){
                if (cpd.dsID!=prevCpd.dsID){
                    Peak pk(cpd.dsID, cpd.cID);
                    pc.push_back(pk);
                }
                prevCpd=cpd;
            }

            if ((int)pc.size()==AppearTimes){
                AlignedCompound ac;
                ac.AlignID=AlignID++;
                ac.Peaks  =pc;
                alignList.push_back(ac);
            }
		}

		alignGroup.clear();
#ifdef _DEBUG
        cout << "press any key to continue the iteration" << endl;
        system("pause");
        cout<<endl<<endl;
#endif
	}

#ifdef _TEST
    ofstream fout;
    fout.open("aligntest",ofstream::out);
    for(auto alignCpd:alignList){
        fout<<alignCpd.toString(true);
    }
    fout.close();
#endif

    ///Output Aligned Result
    /// CSV Title:
    /// AlignID, ds1.mz, ds1.rt,     ds1.sample1.CRR, ds1.sample2.CRR, ds1.sample3.HCC, ds1.sample4.HCC,   ds1.isotopes, ds1.adduct
    ///          ds2.mz, ds2.rt,     ds2.sample1.CRR, ds2.sample2.CRR, ds2.sample3.HCC, ds2.sample4.HCC,   ds2.isotopes, ds2.adduct
    ///          ds3.mz, ds3.rt,     ds3.sample1.CRR, ds3.sample2.CRR, ds3.sample3.HCC, ds3.sample4.HCC,   ds3.isotopes, ds3.adduct
    ///                              format: [ds_name].[samplename].[remark]
    // Title generation
    string title="AlignID,";
    for (size_t dsID=0; dsID!=dsGroups.size();++dsID) {
        size_t datBeginIdx=datBeginIndices[dsID];
        size_t datEndIdx  =datEndIndices[dsID];
        string ds_name="ds"+i2s(dsID);
        title+= ds_name+".mz," + ds_name+".rt,";
        for (size_t col=datBeginIdx; col<=datEndIdx; ++col){
            string samplename = dsGroups[dsID].table(0,col);
            string remark     = dsGroups[dsID].table(1,col);
            title+= ds_name+"."+samplename+"."+remark+",";
        }
        title+= ds_name+".isotopes," + ds_name+".adduct,";
    }

    // Context filling

    ofstream csvOut;
    csvOut.open("align_result.csv", ofstream::out);
    csvOut<<title<<endl;
    for (auto alignCpd:alignList){
        csvOut<<alignCpd.AlignID<<",";
        //cout<<alignCpd.AlignID<<",";
        string datline="";
        for (auto pk: alignCpd.Peaks){
            size_t datBeginIdx=datBeginIndices[pk.dsID];
            size_t datEndIdx  =datEndIndices[pk.dsID];
            for (size_t col=0; col<=datEndIdx+2; ++col){ //copy all columns
                datline+= ("\""+dsGroups[pk.dsID].table(pk.cID+2, col)+"\"" +",");
            }
            //cout<<datline<<endl;
            //char a;cin>>a;
        }
        //datline.erase(datline.size()-1 ,1);
        csvOut<<datline<<endl;

    }
    csvOut.close();


/*
    //Context filling
    ofstream csvOut;
    csvOut.open("align_result.csv",ofstream::out);
    csvOut<<title<<"\n";
    for (auto alignCpd:alignList){
        csvOut<< alignCpd.AlignID<<",";
        for (auto pk:alignCpd.Peaks){
            size_t datBeginIdx= datBeginIndices[pk.dsID];
            size_t datEndIdx  = datEndIndices[pk.dsID];
            for (size_t col=0; col<=datEndIdx+2; ++col)
                csvOut<< "\""+dsGroups[pk.dsID].table(pk.cID+2, col)+"\"" <<",";
        }
        csvOut<<endl;
    }
    csvOut.close();
*/
    return 0;
}

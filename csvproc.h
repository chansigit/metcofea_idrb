#ifndef CSVPROC_H_INCLUDED
#define CSVPROC_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>
#include <cstring>
using namespace std;

class CsvFile{
public:
    ifstream ifs;
    char*    filename;
    vector<vector<string>> grid;
public:
    CsvFile() {}
    CsvFile(char* filename);
    CsvFile(const CsvFile&);
    ~CsvFile();
    string table (size_t x, size_t y);
    vector<string> row   (size_t x);
    vector<string> column(size_t y);
    size_t xsize() {return grid.size();}

    size_t rowcnt(size_t x);
    size_t colcnt(size_t y);

private:
    void read();
    void getCsvLine(string line, vector<string>& token);
};

string CsvFile::table(size_t x, size_t y){
    if (  !( x<grid.size()      && x>=0)  ) return "ERROR";
    if (  !( y<(grid[x]).size() && y>=0)  ) return "ERROR";
    return grid[x][y];
}

CsvFile::CsvFile(char* name){
    this->filename = new char[strlen(name)+1];
    strcpy(filename, name);
    ifs.open(filename, ifstream::in);
    read();
}

CsvFile::CsvFile(const CsvFile& csv){
    this->filename = new char[strlen(csv.filename)+1];
    strcpy(filename, csv.filename);
    ifs.open(filename,ifstream::in);
    read();
}

CsvFile::~CsvFile(){
    ifs.close();
}



void CsvFile::getCsvLine(string line, vector<string>& token){
    string tmp;
    size_t prepos = 0, nextpos = 0;
    while ((nextpos = line.find_first_of(",", prepos)) != string::npos) {
        tmp = line.substr(prepos, nextpos - prepos);
        prepos = nextpos + 1;
        token.push_back(tmp);
    }
    tmp = line.substr(prepos, line.length() - prepos);
    token.push_back(tmp);//for the last section
}

void CsvFile::read(){
    string line;
    while (getline(ifs,line)){
        vector<string> tokens;
        getCsvLine(line,tokens);
        grid.push_back(tokens);
    }
}

vector<string> CsvFile::row   (size_t x){
    return grid[x];
}

vector<string> CsvFile::column(size_t y){
    vector<string> ret;
    for(size_t x=0;x!=grid.size();++x){
        if (y<grid[x].size())
            ret.push_back(grid[x][y]);
        else
            ret.push_back("NULL");
    }
    return ret;
}

size_t CsvFile::rowcnt(size_t x){
    return grid[x].size();
}

size_t CsvFile::colcnt(size_t y){
    size_t cnt=0;
    for(size_t x=0;x!=grid.size();++x){
        if (y<grid[x].size())
            ++cnt;
    }
    return cnt;
}


#endif // CSVPROC_H_INCLUDED

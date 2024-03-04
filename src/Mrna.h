//
// Created by Edward on 2022/4/24.
//

#ifndef HGVSGO_MRNA_H
#define HGVSGO_MRNA_H

#include <unordered_map>
#include <string>

using namespace std;

class Mrna {

public:
    Mrna(const string &mrnaFile);

    string GetSeq(string transcript_id, int begin, int end);
    string GetSeq(string transcript_id, int begin);
    int GetMrnaLength(string transcript_id);

private:
    string mrna_file;
    unordered_map<string, string> transcript_id2seq;
    void ParseMrnaFile(const string &mrnaFile);

};


#endif //HGVSGO_MRNA_H

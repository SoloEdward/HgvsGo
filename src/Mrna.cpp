//
// Created by Edward on 2022/4/24.
//

#include "Mrna.h"
#include "Util.h"
#include <fstream>
#include <iostream>

Mrna::Mrna(const string &mrnaFile) {
    this->ParseMrnaFile(mrnaFile);
}

string Mrna::GetSeq(string transcript_id, int begin, int end) {
    string &seq = this->transcript_id2seq[transcript_id];
    return seq.substr(begin, end - begin);
}

void Mrna::ParseMrnaFile(const string &mrnaFile) {
    ifstream file;
    file.open(mrnaFile);
    if (!file.is_open()) {
        throw runtime_error("Unable to open the genome file");
    }
    string line;
    string seqName;
    string seq;
    long offset = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (seqName.length() > 0) {
                this->transcript_id2seq[seqName] = seq;
            }
            seqName = Util::stringSplit(line.substr(1, line.length()), ' ')[0];
            seq = "";
        } else {
            seq = seq + line;
        }
    }
    this->transcript_id2seq[seqName] = seq;
}

string Mrna::GetSeq(string transcript_id, int begin) {
    string &seq = this->transcript_id2seq[transcript_id];
    return seq.substr(begin);
}

int Mrna::GetMrnaLength(string transcript_id) {
    string &seq = this->transcript_id2seq[transcript_id];
    return seq.size();
}

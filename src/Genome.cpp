//
// Created by Edward on 2022/4/22.
//

#include "Genome.h"
#include <fstream>
#include <iostream>
#include <sstream>

Genome::Genome(const string &genomePath) {
    this->ParseGenomeFile(genomePath);
};


void Genome::ParseGenomeFile(const string &genomeFilePath) {
    ifstream file;
    file.open(genomeFilePath);
    if (!file.is_open()) {
        throw runtime_error("Unable to open the genome file");
    }
    string line;
    string seqName;
    long offset = 0;
    while (getline(file, line)) {
        if (line[0] == '>') {
            seqName = line.substr(1, line.length());
            chrom2offset[seqName] = offset;
        } else {
            for (auto &item : line) {
                genomeArr[offset] = item;
                offset++;
            }
        }
    }
}


char Genome::fetch(string chrom, int pos) { // 1-based
    long offset = this->chrom2offset[chrom];
    return genomeArr[(long) pos - 1 + offset];
}

string Genome::fetch(string chrom, int begin, int end) { // 0-based
    stringstream ss;
    long offset = this->chrom2offset[chrom];
    for (int i = begin; i < end; i++) {
        ss << genomeArr[(long) i + offset];
    }
    return ss.str();
}

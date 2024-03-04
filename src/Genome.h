//
// Created by Edward on 2022/4/22.
//

#ifndef HGVSGO_GENOME_H
#define HGVSGO_GENOME_H

#include <vector>
#include <string>
#include <map>
#include <memory>

using namespace std;

class Genome {

public:
    Genome(const string& genomePath);

    char fetch(string chrom, int pos);  // 1-based pos;
    string fetch(string chrom, int begin, int end);  // 0-based range;

private:
    std::unique_ptr<char[]> genomeArr{new char[4000000000]};
    map<string, long> chrom2offset;

private:
    void ParseGenomeFile(const string &genomeFilePath);
};


#endif //HGVSGO_GENOME_H

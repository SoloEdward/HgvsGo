//
// Created by Edward on 2022/4/25.
//

#ifndef HGVSGO_TRANSCRIPTINTERVALFOREST_H
#define HGVSGO_TRANSCRIPTINTERVALFOREST_H

#include <unordered_map>

#include "Transcript.h"
#include "intervaltree.hpp"

using namespace std;
using namespace Intervals;

class TranscriptIntervalForest {

public:
    TranscriptIntervalForest(string transcriptFile);


    vector<Transcript *> GetTranscripts(string chrom, int begin, int end);

private:
    unordered_map<string, IntervalTree<int, Transcript *>> intervalForest;

    void parseTranscriptFile(string transcriptFile);

};


#endif //HGVSGO_TRANSCRIPTINTERVALFOREST_H

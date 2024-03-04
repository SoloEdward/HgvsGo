//
// Created by Edward on 2022/4/22.
//

#ifndef HGVSGO_TRANSCRIPT_H
#define HGVSGO_TRANSCRIPT_H

#include <string>
#include <iostream>
#include <vector>
#include "intervaltree.hpp"
#include "Exon.h"
#include "Intron.hpp"

using namespace std;

struct CdsPosResult {

    CdsPosResult(const string &cdsPosString, int cdsPos, int exonId, bool isInsideTranscript, bool isOnExon,
                 bool isInsideCds);

    string cdsPosString;
    int cdsPos;
    int exonId;
    bool isInsideTranscript;
    bool isOnExon;
    bool isOnIntron;
    bool isInsideCds;
};

class Transcript {

public:
    string transcript_id;
    string gene;
    bool isReverse;
    string chrom;
    int transcript_start;
    int transcript_end;
    int cds_start;
    int cds_end;
    int cds_length;
    vector<int> cdsOffsets;
    int cdsStartOffset;
    int cdsEndOffset;
    int cdsStartLengthOffset;
    int cdsEndLengthOffset;
    Intervals::IntervalTree<int, Exon *> exons_tree;
    Intervals::IntervalTree<int, Intron *> introns_tree;

public:
    Transcript(const string &transcript_id, string gene, bool isReverse, string chrom, int transcript_start,
               int transcript_end, int cds_start, int cds_end, const string &exon_starts_string,
               const string &exon_ends_string);

    // return cdsPosString,  cdsPos, exonId,
    CdsPosResult GetCdsPosString(int targetPos); // Notice that targetPos is 0-based;

private:
    tuple<string, int, int> GetCdsPosStringOnExon(int targetPos); // return hgvsC, exonId, cdsPos

    tuple<string, int, int> GetCdsPosStringOnIntron(int targetPos);

    tuple<string, int, int> GetCdsPosStringOutsideTranscript(int targetPos);

private:

    vector<int> CalculateExonLengths(vector<int> &exon_starts, vector<int> &exon_ends);

    vector<Exon *> InitExons(vector<int> &exon_starts, vector<int> &exon_ends);

    vector<Intron *> InitIntrons(vector<Exon *> &exons, vector<int> &exon_starts, vector<int> &exon_ends);

    Intervals::IntervalTree<int, Exon *> CreateIntervalTree(vector<Exon *> &exons);

    Intervals::IntervalTree<int, Intron *> CreateIntervalTree(vector<Intron *> &introns);

    int FindExonIndex(Intervals::IntervalTree<int, Exon *> &exons_tree, int targetPos);

    int CalculateCdsLength(vector<int> &exonLengths, int cdsStart, int cdsEnd, vector<Exon *> &exons,
                           Intervals::IntervalTree<int, Exon *> &exons_tree);

    vector<int> CalculateCdsOffsets(vector<int> &exonLengths, int cdsStart, vector<Exon *> &exons,
                                    Intervals::IntervalTree<int, Exon *> &exons_tree);

    int CalculateCdsStartOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsStart);

    int CalculateCdsEndOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsEnd);

    int CalculateCdsStartLengthOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsStart);

    int CalculateCdsEndLengthOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsEnd);

    tuple<int, int> GetCdsPosOnExon(int targetPos); // return cdsPos, exonId

    string FormatCdsString(int cdsPos);
};


#endif //HGVSGO_TRANSCRIPT_H

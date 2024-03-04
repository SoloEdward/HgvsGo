//
// Created by Edward on 2022/4/23.
//

#ifndef HGVSGO_HGVSC_H
#define HGVSGO_HGVSC_H

# include "Variant.h"
# include "Transcript.h"
# include "Genome.h"

struct HgvscResult {

public:

    HgvscResult(const string &type, const string &hgvsC, int exonId, int cdsStartPos, int cdsEndPos,
                const string &altBases, bool isFullyOnIntron, bool isInsideTranscript, bool isInsideCds);

    string type;
    string hgvsC;
    int exonId;
    int cdsStartPos;
    int cdsEndPos;
    string altBases;
    bool isFullyOnIntron;
    bool isInsideTranscript;
    bool isInsideCds;
};

class Hgvsc {

public:
    static HgvscResult AnnotateHgvsc(Variant &v, Transcript &t, Genome &genome);
    // return hgvsC, exonId, cdsStartPos, cdsEndPos  (1-based cdsPos [])

private:
    static HgvscResult AnnotateSingleSnv(Variant &v, Transcript &t);

    static HgvscResult AnnotateDel(Variant &v, Transcript &t);

    static HgvscResult AnnotateIns(Variant &v, Transcript &t, Genome &genome);

    static HgvscResult AnnotateDelIns(Variant &v, Transcript &t);

};


#endif //HGVSGO_HGVSC_H

//
// Created by Edward on 2022/4/24.
//

#ifndef HGVSGO_HGVSP_H
#define HGVSGO_HGVSP_H

#include <string>
#include "Transcript.h"
#include "Translator.h"
#include "Hgvsc.h"
#include "Mrna.h"

using namespace std;

struct HgvspResult {
    HgvspResult(const string &hgvsP);

    string hgvsP;
};

class Hgvsp {

public:
    static HgvspResult
    AnnotateHgvsp(HgvscResult &hgvscResult, Transcript &transcript, Mrna &mrna, Translator &translator);

    static HgvspResult
    ToHgvsp(string &originSeq, string &newSeq, int cdsOffset1, int cdsOffset2, Transcript &transcript, Mrna &mrna,
            Translator &translator, int utrOffset);

    static HgvspResult
    ToFrameShift(string &originSeq, string &newSeq, int cdsOffset1, int cdsOffset2, Transcript &transcript, Mrna &mrna,
                 Translator &translator, int utrOffset);
};


#endif //HGVSGO_HGVSP_H

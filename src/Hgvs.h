//
// Created by Edward on 2022/4/25.
//

#ifndef HGVSGO_HGVS_H
#define HGVSGO_HGVS_H

#include "Variant.h"
#include "Hgvsc.h"
#include "Hgvsp.h"
#include "Transcript.h"
#include "Translator.h"
#include "Mrna.h"
#include "Genome.h"

class Hgvs {
public:
    static tuple<HgvscResult, HgvspResult>
    AnnotateHgvs(Variant &v, Transcript &t, Genome &g, Translator &translator, Mrna &mrna);

};


#endif //HGVSGO_HGVS_H

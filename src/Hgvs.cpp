//
// Created by Edward on 2022/4/25.
//

#include "Hgvs.h"

tuple<HgvscResult, HgvspResult>
Hgvs::AnnotateHgvs(Variant &v, Transcript &t, Genome &g, Translator &translator, Mrna &mrna) {
    HgvscResult hgvscResult = Hgvsc::AnnotateHgvsc(v, t, g);
    HgvspResult hgvspResult = Hgvsp::AnnotateHgvsp(hgvscResult, t, mrna, translator);
    return make_tuple(hgvscResult, hgvspResult);
}

//
// Created by Edward on 2022/4/23.
//

#include "Hgvsc.h"
#include "SeqUtil.h"

HgvscResult::HgvscResult(const string &type, const string &hgvsC, int exonId, int cdsStartPos, int cdsEndPos,
                         const string &altBases, bool isFullyOnIntron, bool isInsideTranscript, bool isInsideCds)
        : type(type), hgvsC(hgvsC), exonId(exonId), cdsStartPos(cdsStartPos), cdsEndPos(cdsEndPos), altBases(altBases),
          isFullyOnIntron(isFullyOnIntron), isInsideTranscript(isInsideTranscript), isInsideCds(isInsideCds) {}

HgvscResult Hgvsc::AnnotateHgvsc(Variant &v, Transcript &t, Genome &g) {
    if (v.ref.length() == 1 and v.alt.length() == 1) {
        return AnnotateSingleSnv(v, t);
    }
    if (v.alt.length() == 1 and v.ref.length() > 1 and v.ref[0] == v.alt[0]) {
        return AnnotateDel(v, t);
    }
    if (v.ref.length() == 1 and v.alt.length() > 1 and v.ref[0] == v.alt[0]) {
        return AnnotateIns(v, t, g);
    }
    return AnnotateDelIns(v, t);
}


HgvscResult Hgvsc::AnnotateSingleSnv(Variant &v, Transcript &t) {
    auto cdsPosResult = t.GetCdsPosString(v.begin);
    string cdsPosString = cdsPosResult.cdsPosString;
    int cdsPos = cdsPosResult.cdsPos;
    int exonId = cdsPosResult.exonId;
    if (t.isReverse) {
        string altBases = SeqUtil::reverseComplement(v.alt);
        return HgvscResult("SNV",
                           "c." + cdsPosString + SeqUtil::reverseComplement(v.ref) + ">" +
                           altBases,
                           exonId, cdsPos, cdsPos, altBases, cdsPosResult.isOnIntron, cdsPosResult.isInsideTranscript,
                           cdsPosResult.isInsideCds);
    }
    return HgvscResult("SNV", "c." + cdsPosString + v.ref + ">" + v.alt, exonId, cdsPos, cdsPos, v.alt,
                       cdsPosResult.isOnIntron, cdsPosResult.isInsideTranscript, cdsPosResult.isInsideCds);
}

HgvscResult Hgvsc::AnnotateDel(Variant &v, Transcript &t) {
    string delBases = v.ref.substr(v.alt.length());

    int delLength = delBases.length();

    auto cdsPosResult1 = t.GetCdsPosString(v.begin + v.alt.length());
    string cdsPosString1 = cdsPosResult1.cdsPosString;
    int exonId1 = cdsPosResult1.exonId;
    int cdsPos1 = cdsPosResult1.cdsPos;

    auto cdsPosResult2 = t.GetCdsPosString(v.end - 1);
    string cdsPosString2 = cdsPosResult2.cdsPosString;
    int exonId2 = cdsPosResult2.exonId;
    int cdsPos2 = cdsPosResult2.cdsPos;


    if (t.isReverse) {
        if (delLength == 1) {
            return HgvscResult("DEL", "c." + cdsPosString1 + "del", exonId1, cdsPos1, cdsPos1, "",
                               cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                               cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                               cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
        }
        return HgvscResult("DEL", "c." + cdsPosString2 + "_" + cdsPosString1 + "del", exonId1, cdsPos2, cdsPos1, "",
                           cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    if (delLength == 1) {
        return HgvscResult("DEL", "c." + cdsPosString1 + "del", exonId2, cdsPos1, cdsPos1, "",
                           cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    return HgvscResult("DEL", "c." + cdsPosString1 + "_" + cdsPosString2 + "del", exonId2, cdsPos1, cdsPos2, "",
                       cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                       cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                       cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
}

HgvscResult Hgvsc::AnnotateIns(Variant &v, Transcript &t, Genome &genome) {
    auto cdsPosResult1 = t.GetCdsPosString(v.begin + v.ref.length() - 1);
    string cdsPosString1 = cdsPosResult1.cdsPosString;
    int exonId1 = cdsPosResult1.exonId;
    int cdsPos1 = cdsPosResult1.cdsPos;
    auto cdsPosResult2 = t.GetCdsPosString(v.begin + v.ref.length());
    string cdsPosString2 = cdsPosResult2.cdsPosString;
    int exonId2 = cdsPosResult2.exonId;
    int cdsPos2 = cdsPosResult2.cdsPos;
    if (t.isReverse) {
        string rightBases = genome.fetch(v.chrom, v.begin + v.ref.length(), v.begin + v.alt.length());
        string altBases = SeqUtil::reverseComplement(v.alt.substr(v.ref.length()));
        if (rightBases == v.alt.substr(v.ref.length())) {
            if (rightBases.length() == 1) {
                return HgvscResult("DUP", "c." + cdsPosString2 + "dup", exonId1, cdsPos2, cdsPos2, altBases,
                                   cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                                   cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                                   cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
            }
            cdsPosResult1 = t.GetCdsPosString(v.begin - 1 + v.alt.length());
            cdsPosString1 = cdsPosResult1.cdsPosString;
            exonId1 = cdsPosResult1.exonId;
            cdsPos1 = cdsPosResult1.cdsPos;
            return HgvscResult("DUP", "c." + cdsPosString1 + "_" + cdsPosString2 + "dup", exonId1, cdsPos1, cdsPos2,
                               altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                               cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                               cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
        }
        return HgvscResult("INS",
                           "c." + cdsPosString2 + "_" + cdsPosString1 + "ins" +
                           altBases,
                           exonId1, cdsPos2, cdsPos1, altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    string leftBases = genome.fetch(v.chrom, v.begin + v.ref.length() - (v.alt.length() - v.ref.length()),
                                    v.begin + v.ref.length());
    string altBases = v.alt.substr(v.ref.length());
    if (leftBases == v.alt.substr(v.ref.length())) {
        if (leftBases.length() == 1) {
            return HgvscResult("DUP", "c." + cdsPosString1 + "dup", exonId2, cdsPos1, cdsPos1, altBases,
                               cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                               cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                               cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
        }
        cdsPosResult2 = t.GetCdsPosString(
                v.begin + v.ref.length() - (v.alt.length() - v.ref.length()));
        string cdsPosString2 = cdsPosResult2.cdsPosString;
        int exonId2 = cdsPosResult2.exonId;
        int cdsPos2 = cdsPosResult2.cdsPos;
        return HgvscResult("DUP", "c." + cdsPosString2 + "_" + cdsPosString1 + "dup", exonId2, cdsPos2, cdsPos1,
                           altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    return HgvscResult("INS", "c." + cdsPosString1 + "_" + cdsPosString2 + "ins" + altBases,
                       exonId2,
                       cdsPos1, cdsPos2, altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                       cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                       cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
}

HgvscResult Hgvsc::AnnotateDelIns(Variant &v, Transcript &t) {
    auto cdsPosResult1 = t.GetCdsPosString(v.begin);
    string cdsPosString1 = cdsPosResult1.cdsPosString;
    int exonId1 = cdsPosResult1.exonId;
    int cdsPos1 = cdsPosResult1.cdsPos;
    auto cdsPosResult2 = t.GetCdsPosString(v.end - 1);
    string cdsPosString2 = cdsPosResult2.cdsPosString;
    int exonId2 = cdsPosResult2.exonId;
    int cdsPos2 = cdsPosResult2.cdsPos;
    if (t.isReverse) {
        string refBases = SeqUtil::reverseComplement(v.ref);
        string altBases = SeqUtil::reverseComplement(v.alt);
        string delinsString = "delins" + altBases;
        if (refBases == SeqUtil::reverseComplement(altBases)) {
            delinsString = "inv";
        }
        if (cdsPosString1 == cdsPosString2) {
            return HgvscResult("DELINS", "c." + cdsPosString2 + delinsString, exonId1,
                               cdsPos2,
                               cdsPos1, altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                               cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                               cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
        }
        return HgvscResult("DELINS",
                           "c." + cdsPosString2 + "_" + cdsPosString1 + delinsString,
                           exonId1, cdsPos2, cdsPos1, altBases, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    string delinsString = "delins" + v.alt;
    if (v.ref == SeqUtil::reverseComplement(v.alt)) {
        delinsString = "inv";
    }
    if (cdsPosString1 == cdsPosString2) {
        return HgvscResult("DELINS", "c." + cdsPosString2 + delinsString, exonId2, cdsPos1, cdsPos2, v.alt,
                           cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                           cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                           cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
    }
    return HgvscResult("DELINS", "c." + cdsPosString1 + "_" + cdsPosString2 + delinsString, exonId2, cdsPos1,
                       cdsPos2, v.alt, cdsPosResult1.isOnIntron && cdsPosResult2.isOnIntron,
                       cdsPosResult1.isInsideTranscript && cdsPosResult2.isInsideTranscript,
                       cdsPosResult1.isInsideCds || cdsPosResult2.isInsideCds);
}



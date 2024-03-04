//
// Created by Edward on 2022/4/24.
//

#include "Hgvsp.h"
#include "Hgvsc.h"
#include <cmath>
#include <algorithm>
#include <regex>

HgvspResult::HgvspResult(const string &hgvsP) : hgvsP(hgvsP) {}

HgvspResult Hgvsp::AnnotateHgvsp(HgvscResult &hgvscResult, Transcript &transcript, Mrna &mrna, Translator &translator) {

    if (hgvscResult.isFullyOnIntron) {
        return HgvspResult("NA");
    }
    if (!hgvscResult.isInsideTranscript or !hgvscResult.isInsideCds) {
        return HgvspResult("NA");
    }
    if (hgvscResult.type == "DUP") {
        std::regex regex1("c\\.[0-9]+_[0-9]+\\+[0-9]+dup");
        bool matched = std::regex_match(hgvscResult.hgvsC, regex1);
        if (matched) {
            return HgvspResult("NA");
        }
    }

    int startCdsPos = hgvscResult.cdsStartPos;
    int endCdsPos = hgvscResult.cdsEndPos;


    int utrOffset = transcript.cdsStartLengthOffset;
    if (transcript.isReverse) {
        utrOffset = transcript.cdsEndLengthOffset;
        startCdsPos = transcript.cds_length + 1 - startCdsPos;
        endCdsPos = transcript.cds_length + 1 - endCdsPos;
    }

    int minValue = min(startCdsPos, endCdsPos);
    int maxValue = max(startCdsPos, endCdsPos);
    if (hgvscResult.type == "INS") {
        startCdsPos = minValue;
    } else {
        startCdsPos = minValue - 1;
    }

    endCdsPos = maxValue;

    if (startCdsPos < 0 or endCdsPos > transcript.cds_length) {
        return HgvspResult("NULL");
    }

    string altBases = hgvscResult.altBases;
    int altLength = altBases.length();
    int mutationLength = 0;


    if (hgvscResult.type == "DELINS") {
        if (altLength > (endCdsPos - startCdsPos)) {
            mutationLength = altBases.length() - (endCdsPos - startCdsPos);
        } else {
            mutationLength = endCdsPos - startCdsPos - altBases.length();
        }
    } else if (hgvscResult.type == "INS" or hgvscResult.type == "DUP") {
        mutationLength = altLength;
    } else if (hgvscResult.type == "DEL") {
        mutationLength = endCdsPos - startCdsPos;
    }

    string hgvspType;
    if (mutationLength % 3 == 0) {
        hgvspType = "INFRAME";
    } else {
        hgvspType = "FRAMESHIFT";
    }
    int cdsOffset1 = floor((double) startCdsPos / 3) * 3;
    int cdsOffset2 = ceil((double) endCdsPos / 3) * 3;

    string originSeq = mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset1,
                                   utrOffset + cdsOffset2);
    string newSeq;

    if (hgvscResult.type == "DUP") {
        newSeq = originSeq.substr(0, endCdsPos - cdsOffset1) + altBases + originSeq.substr(endCdsPos - cdsOffset1);
    } else if (hgvscResult.type == "INS") {
        newSeq = originSeq.substr(0, startCdsPos - cdsOffset1) + altBases + originSeq.substr(startCdsPos - cdsOffset1);
    } else {
        newSeq = originSeq.substr(0, startCdsPos - cdsOffset1) + altBases + originSeq.substr(endCdsPos - cdsOffset1);
    }

    if (hgvspType == "INFRAME") {
        return Hgvsp::ToHgvsp(originSeq, newSeq, cdsOffset1, cdsOffset2, transcript, mrna, translator, utrOffset);
    } else {
        return Hgvsp::ToFrameShift(originSeq, newSeq, cdsOffset1, cdsOffset2, transcript, mrna, translator, utrOffset);
    }
}

HgvspResult
Hgvsp::ToHgvsp(string &originSeq, string &newSeq, int cdsOffset1, int cdsOffset2, Transcript &transcript, Mrna &mrna,
               Translator &translator, int utrOffset) {
    vector<string> originAas = translator.TranslateAas(originSeq);
    vector<string> newAas = translator.TranslateAas(newSeq);


    while (originAas.size() > 0 and newAas.size() > 0 and (originAas[0] == newAas[0]) and
           (originAas.size() > 1 or newAas.size() > 1)) {
        originAas.erase(originAas.begin());
        newAas.erase(newAas.begin());
        cdsOffset1 += 3;
    }

    while (originAas.size() > 0 and newAas.size() > 0 and
           originAas[originAas.size() - 1] == newAas[newAas.size() - 1] and
           (originAas.size() > 1 or newAas.size() > 1)) {
        originAas.erase(originAas.end() - 1);
        newAas.erase(newAas.end() - 1);
    }


    // SNV
    if (originAas.size() == 1 and newAas.size() == 1 and originAas[0] == "Met" and newAas[0] != "MET" and
        cdsOffset1 == 0) {
        return HgvspResult("p.Met1?");
    }

    if (originAas.size() == 1 and newAas.size() == 1 and originAas[0] != "Ter") { // SNV
        if (originAas[0] == newAas[0]) {
            return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "=");
        }
        return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + newAas[0]);
    }

    if (originAas.size() == 1 and newAas.size() == 1 and originAas[0] == "Ter" and newAas[0] != "Ter") { // FRAMESHIFT
        return Hgvsp::ToFrameShift(originSeq, newSeq, cdsOffset1, cdsOffset2, transcript, mrna, translator, utrOffset);
    }

    if (newAas.size() == 0) { // DELETION
        string nextCodon;
        while (cdsOffset1 < transcript.cds_length - 3 * originAas.size()) {
            nextCodon = mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset1 + 3 * originAas.size(),
                                    utrOffset + cdsOffset1 + 3 * originAas.size() + 3);
            string nextAa = translator.TranslateAa(nextCodon);
            if (nextAa == originAas[0]) {
                originAas.erase(originAas.begin());
                originAas.push_back(nextAa);
                cdsOffset1 += 3;
            } else {
                break;
            }
        }
        if (originAas.size() == 1) {
            return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "del");
        }
        return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "_" + originAas[originAas.size() - 1] +
                           to_string(cdsOffset1 / 3 + originAas.size()) + "del");
    }
    if (originAas.size() == 0) { //INSERTION

        string nextCodon;
        int mrnaLength = mrna.GetMrnaLength(transcript.transcript_id);
        while (cdsOffset1 < transcript.cds_length - 3 * newAas.size()) {
            if (utrOffset + cdsOffset1 + 3 * newAas.size() + 3 > mrnaLength){
                break;
            }
            nextCodon = mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset1 + 3 * newAas.size(),
                                    utrOffset + cdsOffset1 + 3 * newAas.size() + 3);
            string nextAa = translator.TranslateAa(nextCodon);
            if (nextAa == newAas[0]) {
                newAas.erase(newAas.begin());
                newAas.push_back(nextAa);
                cdsOffset1 += 3;
            } else {
                break;
            }
        }

        string beforeCodon = mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset1 - 3,
                                         utrOffset + cdsOffset1);
        string beforeAa = translator.TranslateAa(beforeCodon);
        string afterCodon = mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset1,
                                        utrOffset + cdsOffset1 + 3);
        string afterAa = translator.TranslateAa(afterCodon);

        if (newAas.size() == 1) { // only ins 1 aa
            if (newAas[0] == beforeAa) { // 1 aa dup
                return HgvspResult("p." + beforeAa + to_string(cdsOffset1 / 3) + "dup");
            }
            return HgvspResult("p." + beforeAa + to_string(cdsOffset1 / 3) + "_" + afterAa +
                               to_string(cdsOffset1 / 3 + 1) + "ins" + newAas[0]);

        } else {
            int beforeOffset = (int) utrOffset + (int) cdsOffset1 - 3 * (int) newAas.size();
            if (beforeOffset >= 0) {
                string beforeCodons = mrna.GetSeq(transcript.transcript_id,
                                                  utrOffset + cdsOffset1 - 3 * newAas.size(),
                                                  utrOffset + cdsOffset1);
                vector<string> beforeAas = translator.TranslateAas(beforeCodons);
                if (equal(beforeAas.begin(), beforeAas.end(), newAas.begin())) {
                    return HgvspResult("p." + beforeAas[0] + to_string(cdsOffset1 / 3 + 1 - newAas.size()) + "_" +
                                       beforeAas[beforeAas.size() - 1] +
                                       to_string(cdsOffset1 / 3) + "dup");
                }
            }
            string newAasString;
            for (auto &aa : newAas) {
                if (aa == "Ter") {
                    newAasString += aa;
                    break;
                } else {
                    newAasString += aa;
                }
            }
            return HgvspResult("p." + beforeAa + to_string(cdsOffset1 / 3) + "_" +
                               afterAa + to_string(cdsOffset1 / 3 + 1) + "ins" + newAasString);
        }
    }

    if (originAas.size() == 1) {  // DELINS
        if (newAas[0] == "Ter") {
            if (originAas[0] == "Ter") {
                return HgvspResult("p.Ter" + to_string(cdsOffset1 / 3 + 1) + "=");
            }
            return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "Ter");
        }
        string newAasString;
        for (auto &aa : newAas) {
            if (aa == "Ter") {
                newAasString += aa;
                break;
            } else {
                newAasString += aa;
            }
        }
        return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "delins" + newAasString);
    } else {
        string newAasString;
        for (auto &aa : newAas) {
            if (aa == "Ter") {
                newAasString += aa;
                break;
            } else {
                newAasString += aa;
            }
        }
        return HgvspResult("p." + originAas[0] + to_string(cdsOffset1 / 3 + 1) + "_" + originAas[originAas.size() - 1] +
                           to_string(cdsOffset1 / 3 + originAas.size()) + "delins" + newAasString);
    }
}

HgvspResult
Hgvsp::ToFrameShift(string &originSeq, string &newSeq, int cdsOffset1, int cdsOffset2, Transcript &transcript,
                    Mrna &mrna, Translator &translator, int utrOffset) {
    int aaOffset = cdsOffset1 / 3;

    string originCdsSeq = originSeq + mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset2);
    string newCdsSeq = newSeq + mrna.GetSeq(transcript.transcript_id, utrOffset + cdsOffset2);

    int sameAaCount = 0;
    string beforeAa;
    string afterAa;
    for (int i = 0; i < originCdsSeq.size(); i += 3) {
        if (originCdsSeq.substr(i, 3).size() < 3 or newCdsSeq.substr(i, 3).size() < 3) {
            break;
        }
        beforeAa = translator.TranslateAa(originCdsSeq.substr(i, 3));
        afterAa = translator.TranslateAa(newCdsSeq.substr(i, 3));
        if (beforeAa == afterAa and beforeAa != "Ter") {
            sameAaCount++;
        } else {
            break;
        }
    }
    int toTerAaCount = 1;
    string newAa;
    for (int i = 0; i < newCdsSeq.size(); i += 3) {
        if (newCdsSeq.substr(i, 3).size() < 3) {
            break;
        }
        newAa = translator.TranslateAa(newCdsSeq.substr(i, 3));
        if (newAa != "Ter") {
            toTerAaCount++;
        } else {
            break;
        }
    }


    if (beforeAa == "Ter" and afterAa == "Ter") {
        return HgvspResult("p.Ter" + to_string(aaOffset + 1 + sameAaCount) + "=");
    }
    if (beforeAa == "Ter" and afterAa != "Ter") {
        return HgvspResult("p." + beforeAa + to_string(aaOffset + 1 + sameAaCount) + afterAa + "extTer" +
                           to_string(toTerAaCount - sameAaCount - 1));
    }
    if (beforeAa != "Ter" and afterAa == "Ter") {
        return HgvspResult("p." + beforeAa + to_string(aaOffset + 1 + sameAaCount) + "Ter");
    }
    if (newAa != "Ter") { // here we do not find stop codon on the new seq;
        return HgvspResult("p." + beforeAa + to_string(aaOffset + 1 + sameAaCount) + afterAa + "fsTer?");
    }
    return HgvspResult("p." + beforeAa + to_string(aaOffset + 1 + sameAaCount) + afterAa + "fsTer" +
                       to_string(toTerAaCount - sameAaCount));
}




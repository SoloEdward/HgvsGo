//
// Created by Edward on 2022/4/22.
//

#include "Transcript.h"

#include <utility>
#include <cmath>
#include "Util.h"
#include "Exon.h"
#include "Timer.h"

CdsPosResult::CdsPosResult(const string &cdsPosString, int cdsPos, int exonId, bool isInsideTranscript, bool isOnExon,
                           bool isInsideCds)
        : cdsPosString(cdsPosString), cdsPos(cdsPos), exonId(exonId), isInsideTranscript(isInsideTranscript),
          isOnExon(isOnExon), isInsideCds(isInsideCds) {
    this->isOnIntron = isInsideTranscript and not isOnExon;
}

Transcript::Transcript(const string &transcript_id, string gene, bool isReverse, string chrom, int transcript_start,
                       int transcript_end, int cds_start, int cds_end, const string &exon_starts_string,
                       const string &exon_ends_string) {
    this->transcript_id = transcript_id;
    this->gene = std::move(gene);
    this->isReverse = isReverse;
    this->chrom = std::move(chrom);
    this->transcript_id = transcript_id;
    this->transcript_start = transcript_start;
    this->transcript_end = transcript_end;
    this->cds_start = cds_start;
    this->cds_end = cds_end;


    // init exon_starts
    vector<int> exon_starts;
    vector<string> _exon_starts = Util::stringSplit(exon_starts_string, ',');
    for (auto &exon_start: _exon_starts) {
        exon_starts.emplace_back(stoi(exon_start));
    }
    // init exon_ends;
    vector<int> exon_ends;
    vector<string> _exon_ends = Util::stringSplit(exon_ends_string, ',');
    for (auto &exon_end : _exon_ends) {
        exon_ends.emplace_back(stoi(exon_end));
    }

    vector<int> exonLengths = this->CalculateExonLengths(exon_starts, exon_ends);
    vector<Exon *> exons = this->InitExons(exon_starts, exon_ends);
    vector<Intron *> introns = this->InitIntrons(exons, exon_starts, exon_ends);
    this->exons_tree = this->CreateIntervalTree(exons);
    this->introns_tree = this->CreateIntervalTree(introns);
    this->cds_length = this->CalculateCdsLength(exonLengths, this->cds_start, this->cds_end, exons, exons_tree);
    this->cdsOffsets = this->CalculateCdsOffsets(exonLengths, this->cds_start, exons, this->exons_tree);
    this->cdsStartOffset = this->CalculateCdsStartOffset(exon_starts, exon_ends, this->cds_start);
    this->cdsEndOffset = this->CalculateCdsEndOffset(exon_starts, exon_ends, this->cds_end);
    this->cdsStartLengthOffset = this->CalculateCdsStartLengthOffset(exon_starts, exon_ends, this->cds_start);
    this->cdsEndLengthOffset = this->CalculateCdsEndLengthOffset(exon_starts, exon_ends, this->cds_end);
}


vector<int> Transcript::CalculateExonLengths(vector<int> &exon_starts, vector<int> &exon_ends) {
    vector<int> results;
    for (int i = 0; i < exon_starts.size(); i++) {
        results.push_back(exon_ends[i] - exon_starts[i]);
    }
    return results;
}

vector<Exon *> Transcript::InitExons(vector<int> &exon_starts, vector<int> &exon_ends) {
    // init exons;
    vector<Exon *> exons;
    for (int i = 0; i < exon_starts.size(); i++) {
        int s = exon_starts[i];
        int e = exon_ends[i];
        Exon *exon = new Exon(i, s, e);
        exons.push_back(exon);
    }
    return exons;
}

vector<Intron *> Transcript::InitIntrons(vector<Exon *> &exons, vector<int> &exon_starts, vector<int> &exon_ends) {
    //init introns
    vector<Intron *> introns;
    for (int i = 0; i < exon_starts.size() - 1; i++) {
        int s = exon_ends[i];
        int e = exon_starts[i + 1];
        Intron *intron = new Intron(i, s, e);
        introns.push_back(intron);
    }
    // add pointer to introns
    for (int i = 0; i < introns.size(); i++) {
        auto intron = introns[i];
        auto prevExon = exons[i];
        auto nextExon = exons[i + 1];
        intron->SetPrev(prevExon);
        intron->SetNext(nextExon);
    }
    return introns;
}

Intervals::IntervalTree<int, Exon *> Transcript::CreateIntervalTree(vector<Exon *> &exons) {
    // create interval tree
    Intervals::IntervalTree<int, Exon *> exons_tree;
    for (auto &exon: exons) {
        Intervals::Interval<int, Exon *> interval{exon->start, exon->end - 1, exon};
        exons_tree.insert(interval);
    }
    return exons_tree;
}

Intervals::IntervalTree<int, Intron *> Transcript::CreateIntervalTree(vector<Intron *> &introns) {
    Intervals::IntervalTree<int, Intron *> introns_tree;
    for (auto &intron: introns) {
        Intervals::Interval<int, Intron *> interval{intron->start, intron->end - 1, intron};
        introns_tree.insert(interval);
    }
    return introns_tree;
}

int Transcript::FindExonIndex(Intervals::IntervalTree<int, Exon *> &exons_tree, int targetPos) {
    auto overlap_intervals = exons_tree.findIntervalsContainPoint(targetPos);
    auto target_interval = overlap_intervals[0];
    return target_interval.value->index;
}

int Transcript::CalculateCdsLength(vector<int> &exonLengths, int cdsStart, int cdsEnd, vector<Exon *> &exons,
                                   Intervals::IntervalTree<int, Exon *> &exons_tree) {
    int cdsStartIndex = this->FindExonIndex(exons_tree, cdsStart);
    int cdsEndIndex = this->FindExonIndex(exons_tree, cdsEnd - 1);
    if (cdsStartIndex == cdsEndIndex) {
        return cdsEnd - cdsStart;
    }
    int result = 0;
    for (int i = cdsStartIndex + 1; i < cdsEndIndex; i++) {
        result += exonLengths[i];
    }
    result += exons[cdsStartIndex]->end - cdsStart;
    result += cdsEnd - exons[cdsEndIndex]->start;
    return result;
}

vector<int>
Transcript::CalculateCdsOffsets(vector<int> &exonLengths, int cdsStart, vector<Exon *> &exons,
                                Intervals::IntervalTree<int, Exon *> &exons_tree) {
    int cdsStartExonIndex = this->FindExonIndex(exons_tree, cdsStart);
    int to_sub_value = 0;
    for (int i = 0; i < cdsStartExonIndex; i++) {
        to_sub_value += exonLengths[i];
    }
    to_sub_value += cdsStart - exons[cdsStartExonIndex]->start;

    vector<int> cdsPositions;
    for (int i = 0; i < exonLengths.size(); i++) {
        if (i == 0) {
            cdsPositions.push_back(0 - to_sub_value);
        } else {
            cdsPositions.push_back(cdsPositions[i - 1] + exonLengths[i - 1]);
        }
    }
    return cdsPositions;
}

int Transcript::CalculateCdsStartOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsStart) {
    int result = 0;
    for (int i = 0; i < exon_starts.size(); i++) {
        int exonStart = exon_starts[i];
        int exonEnd = exon_ends[i];
        if (cdsStart >= exonEnd) {
            result += exonEnd - exonStart;
        }
        if (cdsStart >= exonStart and cdsStart < exonEnd) {
            result += cdsStart - exonStart;
        }
    }
    return result;
}

int Transcript::CalculateCdsStartLengthOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsStart) {
    int result = 0;
    for (int i = 0; i < exon_starts.size(); i++) {
        int exonStart = exon_starts[i];
        int exonEnd = exon_ends[i];
        if (cdsStart >= exonEnd) {
            result += exonEnd - exonStart;
        }
        if (cdsStart >= exonStart and cdsStart < exonEnd) {
            result += cdsStart - exonStart;
        }
    }
    return result;
}

int Transcript::CalculateCdsEndLengthOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsEnd) {
    int result = 0;
    for (int i = 0; i < exon_starts.size(); i++) {
        int exonStart = exon_starts[i];
        int exonEnd = exon_ends[i];
        if (cdsEnd <= exonStart) {
            result += exonEnd - exonStart;
        }
        if (cdsEnd >= exonStart and cdsEnd < exonEnd) {
            result += exonEnd - cdsEnd;
        }
    }
    return result;
}


int Transcript::CalculateCdsEndOffset(vector<int> &exon_starts, vector<int> &exon_ends, int cdsEnd) {
    auto overlapIntervals = this->exons_tree.findIntervalsContainPoint(cdsEnd - 1);
    return this->cdsOffsets[overlapIntervals[0].value->index] + cdsEnd - 1 - overlapIntervals[0].value->start + 1;
}

string Transcript::FormatCdsString(int cdsPos) {
    if (cdsPos > this->cdsEndOffset) {
        int utrOffset = abs(cdsPos - this->cdsEndOffset);
        if (this->isReverse) {
            return "-" + to_string(utrOffset);
        }
        return "*" + to_string(utrOffset);
    } else if (this->isReverse) {
        if (cdsPos <= 0) {
            return "*" + to_string(abs(cdsPos - 1));
        }
        return to_string(this->cds_length + 1 - cdsPos);
    } else {
        if (cdsPos <= 0) {
            return "-" + to_string(abs(cdsPos - 1));
        }
        return to_string(cdsPos);
    }
}

tuple<int, int> Transcript::GetCdsPosOnExon(int targetPos) {
    auto overlapIntervals = this->exons_tree.findIntervalsContainPoint(targetPos);
    auto targetInterval = overlapIntervals[0];
    auto exon = targetInterval.value;
    int exonId = -1;
    if (this->isReverse) {
        exonId = this->cdsOffsets.size() - exon->index;
    } else {
        exonId = exon->index + 1;
    }
    return make_tuple(this->cdsOffsets[exon->index] + targetPos + 1 - exon->start, exonId);
}

tuple<string, int, int> Transcript::GetCdsPosStringOnExon(int targetPos) {
    tuple<int, int> cdsPosAndExonId = this->GetCdsPosOnExon(targetPos);
    int cdsPos = get<0>(cdsPosAndExonId);
    int exonId = get<1>(cdsPosAndExonId);
    return make_tuple(this->FormatCdsString(cdsPos), cdsPos, exonId);
}

tuple<string, int, int> Transcript::GetCdsPosStringOnIntron(int targetPos) {
    auto overlapIntervals = this->introns_tree.findIntervalsContainPoint(targetPos);
    auto targetInterval = overlapIntervals[0];
    auto intron = targetInterval.value;
    if ((intron->next->start - targetPos) < (targetPos - (intron->prev->end - 1))) {
        int intronOffset = intron->next->start - targetPos;
        int cdsPos = this->cdsOffsets[intron->index + 1] + 1;
        if (this->isReverse) {
            int exonId = this->cdsOffsets.size() - 1 - intron->index;
            return make_tuple(this->FormatCdsString(cdsPos) + "+" + to_string(intronOffset), cdsPos, exonId);
        }
        int exonId = intron->index + 2;
        return make_tuple(this->FormatCdsString(cdsPos) + "-" + to_string(intronOffset), cdsPos, exonId);
    } else {
        int intronOffset = targetPos - intron->prev->end + 1;
        int cdsPos = this->cdsOffsets[intron->index + 1];
        if (isReverse) {
            int exonId = this->cdsOffsets.size() - intron->index;
            return make_tuple(this->FormatCdsString(cdsPos) + "-" + to_string(intronOffset), cdsPos, exonId);
        }
        int exonId = intron->index + 1;
        return make_tuple(this->FormatCdsString(cdsPos) + "+" + to_string(intronOffset), cdsPos, exonId);
    }
}

tuple<string, int, int> Transcript::GetCdsPosStringOutsideTranscript(int targetPos) {
    if (targetPos >= this->transcript_end) {
        tuple<int, int> cdsPosAndExonId = this->GetCdsPosOnExon(this->transcript_end - 1);
        int cdsPos = get<0>(cdsPosAndExonId);
        int distance = targetPos - (this->transcript_end - 1);
        return make_tuple(this->FormatCdsString(cdsPos + distance), cdsPos, -1);
    }
    if (targetPos < this->transcript_start) {
        tuple<int, int> cdsPosAndExonId = this->GetCdsPosOnExon(this->transcript_start);
        int cdsPos = get<0>(cdsPosAndExonId);
        int distance = this->transcript_start - targetPos;
        return make_tuple(this->FormatCdsString(cdsPos + distance), cdsPos, -1);
    }
}

CdsPosResult Transcript::GetCdsPosString(int targetPos) {
    bool isInsideTranscript = false;
    bool isInsideCds = false;
    bool isOnExon = false;
    if (targetPos < this->transcript_start or targetPos >= this->transcript_end) {
        tuple<string, int, int> result = this->GetCdsPosStringOutsideTranscript(targetPos);
        string cdsPosString = get<0>(result);
        int cdsPos = get<1>(result);
        int exonId = get<2>(result);
        return CdsPosResult(cdsPosString, cdsPos, exonId, isInsideTranscript, isOnExon, isInsideCds);
    }
    isInsideTranscript = true;
    if (targetPos >= this->cds_start and targetPos < this->cds_end) {
        isInsideCds = true;
    }
    auto overlapIntervals = this->exons_tree.findIntervalsContainPoint(targetPos);
    if (overlapIntervals.empty()) {
        tuple<string, int, int> result = GetCdsPosStringOnIntron(targetPos);
        string cdsPosString = get<0>(result);
        int cdsPos = get<1>(result);
        int exonId = get<2>(result);
        return CdsPosResult(cdsPosString, cdsPos, exonId, isInsideTranscript, isOnExon, isInsideCds);
    }
    isOnExon = true;
    tuple<string, int, int> result = GetCdsPosStringOnExon(targetPos);
    string cdsPosString = get<0>(result);
    int cdsPos = get<1>(result);
    int exonId = get<2>(result);
    return CdsPosResult(cdsPosString, cdsPos, exonId, isInsideTranscript, isOnExon, isInsideCds);
}




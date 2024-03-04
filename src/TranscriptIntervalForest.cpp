//
// Created by Edward on 2022/4/25.
//
#include <fstream>
#include "TranscriptIntervalForest.h"
#include "Util.h"


TranscriptIntervalForest::TranscriptIntervalForest(string transcriptFile) {
    this->parseTranscriptFile(transcriptFile);
}

vector<Transcript *> TranscriptIntervalForest::GetTranscripts(string chrom, int begin, int end) {
    auto &tree = this->intervalForest[chrom];
    auto overlapIntervals = tree.findOverlappingIntervals({begin, end - 1}, false);
    vector<Transcript *> results;
    for (auto &o : overlapIntervals) {
        results.push_back(o.value);
    }
    return results;
}

void TranscriptIntervalForest::parseTranscriptFile(string transcriptFile) {
    ifstream file;
    file.open(transcriptFile);
    if (!file.is_open()) {
        throw runtime_error("Unable to open the genome file");
    }
    string line;
    getline(file, line);
    vector<string> header = Util::stringSplit(line, '\t');
    int transcript_id_index = distance(header.begin(), find(header.begin(), header.end(), "name"));
    int strand_index = distance(header.begin(), find(header.begin(), header.end(), "strand"));
    int chrom_index = distance(header.begin(), find(header.begin(), header.end(), "chrom"));
    int tx_start_index = distance(header.begin(), find(header.begin(), header.end(), "txStart"));
    int tx_end_index = distance(header.begin(), find(header.begin(), header.end(), "txEnd"));
    int cds_start_index = distance(header.begin(), find(header.begin(), header.end(), "cdsStart"));
    int cds_end_index = distance(header.begin(), find(header.begin(), header.end(), "cdsEnd"));
    int exon_starts_index = distance(header.begin(), find(header.begin(), header.end(), "exonStarts"));
    int exon_ends_index = distance(header.begin(), find(header.begin(), header.end(), "exonEnds"));
    int gene_index = distance(header.begin(), find(header.begin(), header.end(), "name2"));

    while (getline(file, line)) {
        vector<string> value = Util::stringSplit(line, '\t');
        string transcriptId = value[transcript_id_index];
        bool isReverse = value[strand_index] == "-";
        string chrom = value[chrom_index];
        int transcriptStart = stoi(value[tx_start_index]);
        int transcriptEnd = stoi(value[tx_end_index]);
        int cdsStart = stoi(value[cds_start_index]);
        int cdsEnd = stoi(value[cds_end_index]);
        string exonStarts = value[exon_starts_index];
        exonStarts = exonStarts.substr(0, exonStarts.size() - 1);
        string exonEnds = value[exon_ends_index];
        exonEnds = exonEnds.substr(0, exonEnds.size() - 1);
        string gene = value[gene_index];
        Transcript *t = new Transcript(transcriptId, gene, isReverse, chrom, transcriptStart, transcriptEnd, cdsStart,
                                       cdsEnd, exonStarts,
                                       exonEnds);
        if (this->intervalForest.count(chrom) == 0) {
            IntervalTree<int, Transcript *> tree;
            this->intervalForest[chrom] = tree;
        }
        auto &tree = this->intervalForest[chrom];
        if (gene == "TERT") {
            Interval<int, Transcript *> interval{transcriptStart, transcriptEnd - 1 + 500, t};
            tree.insert(interval);
        } else {
            Interval<int, Transcript *> interval{transcriptStart, transcriptEnd - 1, t};
            tree.insert(interval);
        }
    }
}

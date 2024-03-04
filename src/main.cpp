#include <iostream>
#include "Variant.h"
#include "Genome.h"
#include "Timer.h"
#include "SeqUtil.h"
#include "Translator.h"
#include "Mrna.h"
#include "Hgvs.h"
#include "TranscriptIntervalForest.h"
#include <fstream>
#include "Util.h"

using namespace std;

int main(int argc, char *argv[]) {
    {
        if (argc != 6) {
            cout << "Usage: ./HgvsGo RnaFastaFile humanGenomeFile transcriptFile inputFile outputFile" << endl;
            exit(1);
        }

        string mrnaFile = argv[1];
        string genomeFile = argv[2];
        string transcriptFile = argv[3];
        string inputFile = argv[4];
        string outputFile = argv[5];

        Timer timer;
        Mrna mrna(mrnaFile);
        Genome genome(genomeFile);
        TranscriptIntervalForest forest(transcriptFile);
        Translator translator;

        ofstream outf;
        outf.open(outputFile, ios::trunc);
        if (!outf) {
            cerr << "open file error " << outputFile << endl;
        }

        ifstream file;
        file.open(inputFile);
        if (!file.is_open()) {
            throw runtime_error("Unable to open the input file");
        }
        string line;
        getline(file, line);
        outf << line << "\ttranscript_id\tgene\texon_id\thgvs_c\thgvs_p\n";
        vector<string> header = Util::stringSplit(line, '\t');
        int chrom_index = distance(header.begin(), find(header.begin(), header.end(), "chrom"));
        int pos_index = distance(header.begin(), find(header.begin(), header.end(), "pos"));
        int ref_index = distance(header.begin(), find(header.begin(), header.end(), "ref"));
        int alt_index = distance(header.begin(), find(header.begin(), header.end(), "alt"));
        while (getline(file, line)) {
            vector<string> value = Util::stringSplit(line, '\t');
            string chrom = value[chrom_index];
            int pos = stoi(value[pos_index]);
            string ref = value[ref_index];
            string alt = value[alt_index];
            Variant v{chrom, pos, ref, alt};
            Variant v2 = Variant::LeftTrim(v);
            Variant v3 = Variant::RightTrim(v2);
            Variant vPer5 = Variant::Per5Align(v3, genome);
            Variant vPer3 = Variant::Per3Align(v3, genome);
            auto transcriptPtrs = forest.GetTranscripts(v.chrom, v.begin, v.end);
            if (transcriptPtrs.empty()) {
                outf << line << "\tNA\tNA\tNA\tNA\tNA\n";
            }
            for (auto t: transcriptPtrs) {
                Variant toAnnotateVariant = v;
                if (t->isReverse) {
                    toAnnotateVariant = vPer5;
                } else {
                    toAnnotateVariant = vPer3;
                }
                auto hgvsResult = Hgvs::AnnotateHgvs(toAnnotateVariant, (*t), genome, translator, mrna);
                auto hgvscResult = get<0>(hgvsResult);
                auto hgvspResult = get<1>(hgvsResult);
                outf << line << "\t" << t->transcript_id << "\t" << t->gene << "\t" << hgvscResult.exonId << "\t"
                     << hgvscResult.hgvsC << "\t" << hgvspResult.hgvsP << "\n";
            }
        }
        outf.close();

    }
}

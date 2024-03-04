//
// Created by Edward on 2022/4/24.
//

#include "Translator.h"
#include "iostream"

Translator::Translator() {
    this->codon2aa.insert({"TTT", "Phe"});
    this->codon2aa.insert({"TTT", "Phe"});
    this->codon2aa.insert({"TTC", "Phe"});
    this->codon2aa.insert({"TTA", "Leu"});
    this->codon2aa.insert({"TTG", "Leu"});
    this->codon2aa.insert({"TCT", "Ser"});
    this->codon2aa.insert({"TCC", "Ser"});
    this->codon2aa.insert({"TCA", "Ser"});
    this->codon2aa.insert({"TCG", "Ser"});
    this->codon2aa.insert({"TAT", "Tyr"});
    this->codon2aa.insert({"TAC", "Tyr"});
    this->codon2aa.insert({"TAA", "Ter"});
    this->codon2aa.insert({"TAG", "Ter"});
    this->codon2aa.insert({"TGT", "Cys"});
    this->codon2aa.insert({"TGC", "Cys"});
    this->codon2aa.insert({"TGA", "Ter"});
    this->codon2aa.insert({"TGG", "Trp"});
    this->codon2aa.insert({"CTT", "Leu"});
    this->codon2aa.insert({"CTC", "Leu"});
    this->codon2aa.insert({"CTA", "Leu"});
    this->codon2aa.insert({"CTG", "Leu"});
    this->codon2aa.insert({"CCT", "Pro"});
    this->codon2aa.insert({"CCC", "Pro"});
    this->codon2aa.insert({"CCA", "Pro"});
    this->codon2aa.insert({"CCG", "Pro"});
    this->codon2aa.insert({"CAT", "His"});
    this->codon2aa.insert({"CAC", "His"});
    this->codon2aa.insert({"CAA", "Gln"});
    this->codon2aa.insert({"CAG", "Gln"});
    this->codon2aa.insert({"CGT", "Arg"});
    this->codon2aa.insert({"CGC", "Arg"});
    this->codon2aa.insert({"CGA", "Arg"});
    this->codon2aa.insert({"CGG", "Arg"});
    this->codon2aa.insert({"ATT", "Ile"});
    this->codon2aa.insert({"ATC", "Ile"});
    this->codon2aa.insert({"ATA", "Ile"});
    this->codon2aa.insert({"ATG", "Met"});
    this->codon2aa.insert({"ACT", "Thr"});
    this->codon2aa.insert({"ACC", "Thr"});
    this->codon2aa.insert({"ACA", "Thr"});
    this->codon2aa.insert({"ACG", "Thr"});
    this->codon2aa.insert({"AAT", "Asn"});
    this->codon2aa.insert({"AAC", "Asn"});
    this->codon2aa.insert({"AAA", "Lys"});
    this->codon2aa.insert({"AAG", "Lys"});
    this->codon2aa.insert({"AGT", "Ser"});
    this->codon2aa.insert({"AGC", "Ser"});
    this->codon2aa.insert({"AGA", "Arg"});
    this->codon2aa.insert({"AGG", "Arg"});
    this->codon2aa.insert({"GTT", "Val"});
    this->codon2aa.insert({"GTC", "Val"});
    this->codon2aa.insert({"GTA", "Val"});
    this->codon2aa.insert({"GTG", "Val"});
    this->codon2aa.insert({"GCT", "Ala"});
    this->codon2aa.insert({"GCC", "Ala"});
    this->codon2aa.insert({"GCA", "Ala"});
    this->codon2aa.insert({"GCG", "Ala"});
    this->codon2aa.insert({"GAT", "Asp"});
    this->codon2aa.insert({"GAC", "Asp"});
    this->codon2aa.insert({"GAA", "Glu"});
    this->codon2aa.insert({"GAG", "Glu"});
    this->codon2aa.insert({"GGT", "Gly"});
    this->codon2aa.insert({"GGC", "Gly"});
    this->codon2aa.insert({"GGA", "Gly"});
    this->codon2aa.insert({"GGG", "Gly"});
}

string Translator::TranslateAa(string &codon) {
    return this->codon2aa[codon];
}

vector<string> Translator::TranslateAas(string &codons) {
    vector<string> results;
    for (int i = 0; i < codons.size(); i += 3) {
        string codon = codons.substr(i, 3);
        results.push_back(this->TranslateAa(codon));
    }
    return results;
}

string Translator::TranslateAa(string &&codon) {
    return this->codon2aa[codon];
}


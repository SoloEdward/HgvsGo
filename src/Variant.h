//
// Created by Edward on 2022/4/22.
//

#ifndef HGVSGO_VARIANT_H
#define HGVSGO_VARIANT_H

#include <string>
#include <vector>
#include "Genome.h"

using namespace std;

class Variant {

public:
    Variant(string chrom, int pos, string ref, string alt);

    static Variant LeftTrim(Variant v);

    static Variant RightTrim(Variant v);

    static Variant Per5Align(Variant v, Genome &genome);

    static Variant Per3Align(Variant v, Genome &genome);

    static vector<Variant> CreateVariants(string chrom, int pos, string ref, string alt, Genome &genome);

public:
    string chrom;
    int pos;
    string ref;
    string alt;
    int begin;
    int end;
};




//ostream &operator<<(ostream &out, Variant &v) {
//    out << v.chrom + ":" + to_string(v.pos) + " " + v.ref + ">" + v.alt;
//    return out;
//}

#endif //HGVSGO_VARIANT_H

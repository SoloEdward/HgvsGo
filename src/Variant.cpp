//
// Created by Edward on 2022/4/22.
//

#include "Variant.h"
#include <iostream>

Variant::Variant(string chrom, int pos, string ref, string alt) {
    this->chrom = chrom;
    if (chrom.substr(0, 3) != "chr") {
        this->chrom = "chr" + chrom;
    }
    this->pos = pos;
    this->ref = ref;
    this->alt = alt;
    this->begin = pos - 1;
    this->end = pos - 1 + ref.length();
}


Variant Variant::Per5Align(Variant v, Genome &genome) {
    if (v.alt.length() == 1 and v.ref.length() > 1 and v.alt[0] == v.ref[0]) {
        if (v.ref[0] == v.ref[v.ref.length() - 1]) {
            char left_base = genome.fetch(v.chrom, v.pos - 1);
            string newAlt(1, left_base);
            return Per5Align(
                    Variant(v.chrom, v.pos - 1, newAlt + v.ref.substr(0, v.ref.length() - 1),
                            newAlt),
                    genome);
        }
        return v;
    }
    if (v.ref.length() == 1 and v.alt.length() > 1 and v.alt[0] == v.ref[0]) {
        if (v.alt[v.alt.length() - 1] == v.alt[0]) {
            char left_base = genome.fetch(v.chrom, v.pos - 1);
            string newRef(1, left_base);
            return Per5Align(
                    Variant(v.chrom, v.pos - 1, newRef,
                            newRef + v.alt.substr(0, v.alt.length() - 1)),
                    genome);
        }
        return v;
    }
    return v;
};

Variant Variant::Per3Align(Variant v, Genome &genome) {
    if (v.alt.length() == 1 and v.ref.length() > 1 and v.alt[0] == v.ref[0]) {
        char right_base = genome.fetch(v.chrom, v.pos + v.ref.length());
        if (v.ref[1] == right_base) {
            string newAlt(1, v.ref[1]);
            return Per3Align(
                    Variant(v.chrom, v.pos + 1, v.ref.substr(1, v.ref.length()) + newAlt,
                            newAlt),
                    genome);
        }
        return v;
    }
    if (v.ref.length() == 1 and v.alt.length() > 1 and v.alt[0] == v.ref[0]) {
        char right_base = genome.fetch(v.chrom, v.pos + 1);
        if (v.alt[1] == right_base) {
            string newRef(1, v.alt[1]);
            return Per3Align(
                    Variant(v.chrom, v.pos + 1, newRef,
                            v.alt.substr(1, v.alt.length()) + newRef),
                    genome);
        }
        return v;
    }
    return v;
};

vector<Variant> Variant::CreateVariants(string chrom, int pos, string ref, string alt, Genome &genome) {
    vector<Variant> variants;
    Variant v(chrom, pos, ref, alt);
    Variant v2 = Variant::Per5Align(v, genome);
    Variant v3 = Variant::Per3Align(v, genome);
    variants.push_back(v);
    variants.push_back(v2);
    variants.push_back(v3);
    return variants;
}

Variant Variant::LeftTrim(Variant v) {
    string ref = v.ref;
    string alt = v.alt;
    int pos = v.pos;
    while (ref.length() > 1 and alt.length() > 1 and ref[0] == alt[0]) {
        ref = ref.substr(1);
        alt = alt.substr(1);
        pos = pos + 1;
    }
    return Variant(v.chrom, pos, ref, alt);
}

Variant Variant::RightTrim(Variant v) {
    string ref = v.ref;
    string alt = v.alt;
    while (ref.length() > 1 and alt.length() > 1 and ref[ref.length() - 1] == alt[alt.length() - 1]) {
        ref = ref.substr(0, ref.length() - 1);
        alt = alt.substr(0, alt.length() - 1);
    }
    return Variant(v.chrom, v.pos, v.ref, v.alt);
};

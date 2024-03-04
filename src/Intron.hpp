//
// Created by Edward on 2022/4/23.
//

#ifndef HGVSGO_INTRON_H
#define HGVSGO_INTRON_H


#include "Exon.h"

class Intron {
public:
    int index;
    int start;
    int end;
    Exon *prev;
    Exon *next;

public:
    Intron(int index, int start, int end) {
        this->index = index;
        this->start = start;
        this->end = end;
    };

    void SetPrev(Exon *exonPtr) {
        this->prev = exonPtr;
    };

    void SetNext(Exon *exonPtr) {
        this->next = exonPtr;
    }
};

#endif //HGVSGO_INTRON_H

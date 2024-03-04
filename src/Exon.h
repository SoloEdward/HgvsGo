//
// Created by Edward on 2022/4/23.
//

#ifndef HGVSGO_EXON_H
#define HGVSGO_EXON_H


class Exon {
public:
    int index;
    int start;
    int end;

public:
    Exon(int index, int start, int end) {
        this->index = index;
        this->start = start;
        this->end = end;
    }
};


#endif //HGVSGO_EXON_H

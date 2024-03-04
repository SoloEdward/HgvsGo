//
// Created by Edward on 2022/4/23.
//

#include "SeqUtil.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
using namespace std;

string SeqUtil::reverseComplement(string seq) {
    stringstream ss;
    for (int i = seq.size() - 1; i >= 0; i--) {
        switch (seq[i]) {
            case 'A':
                ss << "T";
                break;
            case 'T':
                ss << 'A';
                break;
            case 'C':
                ss << 'G';
                break;
            case 'G':
                ss << 'C';
                break;
            case 'N':
                ss << 'N';
                break;
            default:
                ss << '?';
                break;
        }
    }
    return ss.str();
}

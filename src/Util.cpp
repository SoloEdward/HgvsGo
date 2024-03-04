//
// Created by Edward on 2022/4/23.
//

#include "Util.h"
#include <sstream>

vector<string> Util::stringSplit(const string &str, char delim) {
    vector<string> results;
    string tmp;
    stringstream ss(str);
    while (getline(ss, tmp, delim)) {
        results.push_back(tmp);
    }
    return results;
};
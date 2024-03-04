//
// Created by Edward on 2022/4/24.
//

#ifndef HGVSGO_TRANSLATOR_H
#define HGVSGO_TRANSLATOR_H

#include <unordered_map>
#include <string>
#include <map>
#include <vector>

using namespace std;



class Translator {

public:
    Translator();

    string TranslateAa(string &&codon);
    string TranslateAa(string &codon);

    vector<string> TranslateAas(string &codons);

private:
    std::unordered_map<std::string, std::string> codon2aa;
};


#endif //HGVSGO_TRANSLATOR_H

//
// Created by karel on 25.12.17.
//

#ifndef EWALDIS_CPP_NOMINALENCODER_H
#define EWALDIS_CPP_NOMINALENCODER_H


#include <set>
#include <map>
#include <string>

class NominalEncoder {
private:
    std::map<std::string, int> value_to_encoding;
    std::map<int, std::string> encoding_to_value;
public:
    double get_encoding(std::string value);
    std::string get_value(double encoding);
};


#endif //EWALDIS_CPP_NOMINALENCODER_H

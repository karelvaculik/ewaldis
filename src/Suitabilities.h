//
// Created by karel on 2.1.18.
//

#ifndef EWALDIS_CPP_SUITABILITIES_H
#define EWALDIS_CPP_SUITABILITIES_H

#include <map>


class Suitabilities {

    const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_positive;
    const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_negative;

public:
    Suitabilities(std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_positive,
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_negative);


    const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & get_edge_pairs_dictionary_positive();
    const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & get_edge_pairs_dictionary_negative();
};


#endif //EWALDIS_CPP_SUITABILITIES_H

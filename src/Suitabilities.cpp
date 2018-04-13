//
// Created by karel on 2.1.18.
//

#include <iostream>
#include "Suitabilities.h"


Suitabilities::Suitabilities(
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_positive,
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_negative)
        : edge_pairs_dictionary_positive(edge_pairs_dictionary_positive),
          edge_pairs_dictionary_negative(edge_pairs_dictionary_negative) {

}


const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & Suitabilities::get_edge_pairs_dictionary_positive()
{
    return edge_pairs_dictionary_positive;
};
const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & Suitabilities::get_edge_pairs_dictionary_negative()
{
    return edge_pairs_dictionary_negative;
};

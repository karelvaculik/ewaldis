//
// Created by karel on 14.1.18.
//

#ifndef EWALDIS_CPP_PATTERNCHECKER_H
#define EWALDIS_CPP_PATTERNCHECKER_H


#include "commonutils.h"
#include "RandomGenerator.h"
#include "DynamicGraph.h"
#include "Pattern.h"

class PatternChecker {
private:
    bool use_vertex_attributes;
    timestamp_t time_unit_secondary;
    int random_walks;
    RandomGenerator * random_generator;
public:
    PatternChecker(bool use_vertex_attributes,
                   timestamp_t time_unit_secondary,
                   int random_walks,
                   RandomGenerator * random_generator);
    double check_pattern_in_instance(DynamicGraph * graph, Pattern * pattern,
                                     std::vector<int> & graph_starting_vertices,
                                     timestamp_t graph_starting_timestamp,
                                     std::vector<double> & edge_weights);
};


#endif //EWALDIS_CPP_PATTERNCHECKER_H

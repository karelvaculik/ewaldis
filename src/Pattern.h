//
// Created by karbal on 29.12.17.
//

#ifndef EWALDIS_CPP_PATTERN_H
#define EWALDIS_CPP_PATTERN_H


#include "PatternEdge.h"
#include "NominalEncoder.h"

class Pattern {
private:
    std::vector<std::pair<int, int>> pattern_vertex_pairs;
    std::vector<bool> directions;
    std::vector<std::vector<timestamp_t>> timestamps;
    std::vector<std::vector<std::map<std::string, double>>> attributes;
    std::vector<double> scores;

public:
    Pattern(std::vector<std::pair<int, int>> pattern_vertex_pairs,
            std::vector<bool> directions,
            std::vector<std::vector<timestamp_t>> timestamps,
            std::vector<std::vector<std::map<std::string, double>>> attributes,
            std::vector<double> scores);

    const std::vector<std::pair<int, int>> & get_pattern_vertex_pairs();
    const std::vector<bool> & get_directions();
    const std::vector<std::vector<timestamp_t>> & get_timestamps();
    const std::vector<std::vector<std::map<std::string, double>>> & get_attributes();
    std::vector<std::vector<std::map<std::string, std::string>>> get_decoded_attributes(std::map<std::string, AttributeType> edge_schema,
                                                                                        NominalEncoder * ne);
    const std::vector<double> & get_scores();
    void clean_empty_instances();

};


#endif //EWALDIS_CPP_PATTERN_H

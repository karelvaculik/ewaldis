//
// Created by karbal on 2.1.18.
//

#include "PatternEdge.h"


PatternEdge::PatternEdge(std::vector<Edge*> pattern_edge, double pattern_score) :
        pattern_edge(pattern_edge),
        pattern_score(pattern_score)
{
}


const std::vector<Edge*> & PatternEdge::get_pattern_edge()
{
    return pattern_edge;
}

double PatternEdge::get_pattern_score()
{
    return pattern_score;
}

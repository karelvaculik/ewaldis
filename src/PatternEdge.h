//
// Created by karbal on 2.1.18.
//

#ifndef EWALDIS_CPP_PATTERNEDGE_H
#define EWALDIS_CPP_PATTERNEDGE_H

#include <vector>
#include "Edge.h"

class PatternEdge
{
private:
    std::vector<Edge*> pattern_edge;
    double pattern_score;

public:
    PatternEdge(std::vector<Edge*> pattern_edge, double pattern_score);

    const std::vector<Edge*> & get_pattern_edge();
    double get_pattern_score();

};


#endif //EWALDIS_CPP_PATTERNEDGE_H

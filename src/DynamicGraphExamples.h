//
// Created by karel on 25.12.17.
//

#ifndef EWALDIS_CPP_DYNAMICGRAPHEXAMPLES_H
#define EWALDIS_CPP_DYNAMICGRAPHEXAMPLES_H


#include <vector>
#include "DynamicGraph.h"
#include "NominalEncoder.h"

class DynamicGraphExamples {
private:
    static Vertex ver(int v_id);
    static Edge edg(int e_id, int f_id, int t_id, std::string label, NominalEncoder * ne, timestamp_t timestamp = 0);
public:
    static DynamicGraph prepareExample1(NominalEncoder * ne);
    static DynamicGraph prepareExample2(NominalEncoder * ne);
};


#endif //EWALDIS_CPP_DYNAMICGRAPHEXAMPLES_H

//
// Created by karel on 12.12.17.
//

#ifndef EWALDIS_CPP_EDGE_H
#define EWALDIS_CPP_EDGE_H

#include <map>
#include <string>
#include "commonutils.h"

class Edge {
private:
    int edge_id;
    int original_edge_id;
    int from_vertex_id;
    int to_vertex_id;
    timestamp_t timestamp;
    bool direction;
    std::map<std::string, double> attributes;

public:
    Edge(int edge_id, int from_vertex_id, int to_vertex_id, timestamp_t timestamp,
         std::map<std::string, double> attributes, bool direction=true,
         int original_edge_id=constants::INVALID_ID);

    Edge create_opposite_edge(bool undirected, int new_edge_id);
    int get_edge_id();
    int get_from_vertex_id();
    int get_to_vertex_id();
    timestamp_t get_timestamp();
    const std::map<std::string, double> & get_attributes();
    bool get_direction();
    int get_original_edge_id();
};


#endif //EWALDIS_CPP_EDGE_H

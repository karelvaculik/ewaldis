//
// Created by karel on 12.12.17.
//

#include <iostream>
#include "Edge.h"

Edge::Edge(int edge_id, int from_vertex_id, int to_vertex_id, timestamp_t timestamp,
           std::map<std::string, double> attributes,
           bool direction, int original_edge_id)
{
    this->edge_id = edge_id;
    this->from_vertex_id = from_vertex_id;
    this->to_vertex_id = to_vertex_id;
    this->timestamp = timestamp;
    this->attributes = attributes;
    this->direction = direction;
    if (original_edge_id != constants::INVALID_ID) {
        this->original_edge_id = original_edge_id;
    } else {
        this->original_edge_id = edge_id;
    }

}

Edge Edge::create_opposite_edge(bool undirected, int new_edge_id)
{
    return Edge(new_edge_id, this->to_vertex_id, this->from_vertex_id,
                    this->timestamp, this->attributes, undirected || !this->direction,
                    this->original_edge_id);
}

int Edge::get_edge_id()
{
    return this->edge_id;
}

int Edge::get_from_vertex_id()
{
    return this->from_vertex_id;
}

int Edge::get_to_vertex_id()
{
    return this->to_vertex_id;
}

timestamp_t Edge::get_timestamp()
{
    return this->timestamp;
}

const std::map<std::string, double> & Edge::get_attributes()
{
    return this->attributes;
}

bool Edge::get_direction()
{
    return this->direction;
}

int Edge::get_original_edge_id()
{
    return this->original_edge_id;
}

//
// Created by karel on 12.12.17.
//

#include "Vertex.h"


Vertex::Vertex(int vertex_id, std::map<std::string, double> attributes)
{
    this->vertex_id = vertex_id;
    this->attributes = attributes;
}

int Vertex::get_vertex_id()
{
    return this->vertex_id;
}

const std::map<std::string, double> & Vertex::get_attributes()
{
    return this->attributes;
}
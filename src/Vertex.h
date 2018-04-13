//
// Created by karel on 12.12.17.
//

#ifndef EWALDIS_CPP_VERTEX_H
#define EWALDIS_CPP_VERTEX_H

#include <map>
#include <string>

class Vertex
{

private:
    int vertex_id;
    std::map<std::string, double> attributes;

public:
    Vertex(int vertex_id, std::map<std::string, double> attributes);
    int get_vertex_id();
    const std::map<std::string, double> & get_attributes();
};


#endif //EWALDIS_CPP_VERTEX_H

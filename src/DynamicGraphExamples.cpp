//
// Created by karel on 25.12.17.
//

#include <iostream>
#include "DynamicGraphExamples.h"

Vertex DynamicGraphExamples::ver(int v_id)
{
    std::map<std::string, double> attributes;
    attributes["label"] = 0;
    return Vertex(v_id, attributes);
}

Edge DynamicGraphExamples::edg(int e_id, int f_id, int t_id, std::string label, NominalEncoder * ne,
    timestamp_t timestamp)
{
    std::map<std::string, double> attributes;
    attributes["label"] = ne->get_encoding(label);
    return Edge(e_id, f_id, t_id, timestamp, attributes);
}

DynamicGraph DynamicGraphExamples::prepareExample1(NominalEncoder * ne)
{
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (int i = 0; i < 12; ++i) {
        vertices.push_back(ver(i));
    }
    edges.push_back(edg(1, 0, 1, "a", ne));
    edges.push_back(edg(2, 0, 2, "d", ne));
    edges.push_back(edg(3, 3, 4, "a", ne));
    edges.push_back(edg(4, 3, 5, "d", ne));
    edges.push_back(edg(5, 6, 7, "a", ne));
    edges.push_back(edg(6, 6, 8, "x", ne));
    edges.push_back(edg(7, 9, 10, "a", ne));
    edges.push_back(edg(8, 9, 11, "x", ne));

    std::map<std::string, AttributeType> vertex_schema;
    std::map<std::string, AttributeType> edge_schema;

    vertex_schema["label"] = AttributeType::NOMINAL;
    edge_schema["label"] = AttributeType::NOMINAL;

    DynamicGraph graph = DynamicGraph(vertices, edges, vertex_schema, edge_schema, true);

    return graph;
}


DynamicGraph DynamicGraphExamples::prepareExample2(NominalEncoder * ne)
{
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    for (int i = 0; i < 12; ++i) {
        vertices.push_back(ver(i));
    }
    edges.push_back(edg(1, 0, 1, "a", ne));
    edges.push_back(edg(2, 0, 2, "d", ne));
    edges.push_back(edg(3, 0, 3, "o", ne, 1.0));
    edges.push_back(edg(4, 4, 5, "a", ne));
    edges.push_back(edg(5, 4, 6, "d", ne));
    edges.push_back(edg(6, 4, 7, "o", ne, 1.0));
    edges.push_back(edg(7, 8, 9, "a", ne));
    edges.push_back(edg(8, 8, 10, "x", ne));
    edges.push_back(edg(9, 11, 12, "a", ne));
    edges.push_back(edg(10, 11, 13, "x", ne));

    std::map<std::string, AttributeType> vertex_schema;
    std::map<std::string, AttributeType> edge_schema;

    vertex_schema["label"] = AttributeType::NOMINAL;
    edge_schema["label"] = AttributeType::NOMINAL;

    DynamicGraph graph = DynamicGraph(vertices, edges, vertex_schema, edge_schema, true);

    return graph;
}
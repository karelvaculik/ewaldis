//
// Created by karel on 25.12.17.
//

#ifndef EWALDIS_CPP_DYNAMICGRAPH_H
#define EWALDIS_CPP_DYNAMICGRAPH_H


#include <vector>
#include "Vertex.h"
#include "Edge.h"
#include "commonutils.h"
#include "NominalEncoder.h"
#include "RandomGenerator.h"

class DynamicGraph {
private:
    std::vector<Vertex> original_vertices;
    std::vector<Edge> original_edges;
    std::map<int, Vertex *> vertices;
    std::map<int, Edge *> edges;
    std::map<std::string, AttributeType> vertex_schema;
    std::map<std::string, AttributeType> edge_schema;
    bool undirected;
    std::map<int, std::vector<Edge *>> adjacency_list;

    void create_adjacency_list();
    DynamicGraph create_subgraph_with_bounded_edges(timestamp_t upper_time_limit, timestamp_t lower_time_limit);

public:
    DynamicGraph(std::vector<Vertex> & vertices, std::vector<Edge> & edges,
                 std::map<std::string, AttributeType> & vertex_schema,
                 std::map<std::string, AttributeType> & edge_schema,
                 bool undirected);
    DynamicGraph(std::map<int, Vertex *> vertices, std::map<int, Edge *> edges,
                 std::map<std::string, AttributeType> vertex_schema,
                 std::map<std::string, AttributeType> edge_schema,
                 bool undirected, std::map<int, std::vector<Edge *>> adjacency_list);

    const std::map<int, Vertex *> & get_vertices();
    const std::map<int, Edge *> & get_edges();
    const std::map<std::string, AttributeType> & get_vertex_schema();
    const std::map<std::string, AttributeType> & get_edge_schema();
    bool is_undirected();
    const std::map<int, std::vector<Edge *>> & get_adjacency_list();
    std::vector<DynamicGraph> create_subgraph_instances(std::vector<timestamp_t> positive_event_times,
                                                        std::vector<timestamp_t> negative_event_times,
                                                        timestamp_t primary_time_unit);

    void identify_edge_events(std::string attr_name, int label, int sample_size, RandomGenerator * random_generator,
                              std::vector<int> & edge_ids, std::vector<std::vector<int>> & vertex_ids,
                              std::vector<timestamp_t> & vertex_timestamps);
};





#endif //EWALDIS_CPP_DYNAMICGRAPH_H

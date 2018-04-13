//
// Created by karel on 25.12.17.
//

#include <iostream>
#include "DynamicGraph.h"


#include "commonutilstemplated.h"


DynamicGraph::DynamicGraph(std::vector<Vertex> & vertices, std::vector<Edge> & edges,
                           std::map<std::string, AttributeType> & vertex_schema,
                           std::map<std::string, AttributeType> & edge_schema,
                           bool undirected)
{
    this->vertex_schema = vertex_schema;
    this->edge_schema = edge_schema;
    this->undirected = undirected;

    // prepare the original vertices and edges (we keep the objects in these vectors)
    for (Vertex &value: vertices) {
        this->original_vertices.push_back(value);
    }
    int m = -1;
    for (Edge &value: edges) {
        this->original_edges.push_back(value);
        if (value.get_edge_id() > m) {
            m = value.get_edge_id();
        }
    }
    m += 1;
    for (Edge &value: edges) {
        Edge opposite_edge = value.create_opposite_edge(undirected, value.get_edge_id() + m);
        this->original_edges.push_back(opposite_edge);
    }

    // prepare maps with pointer to vertices and edges (these are then used in adjacency matrices and subgraphs)
    for (Vertex &value: original_vertices)
    {
        this->vertices[value.get_vertex_id()] = &value;
    }
    for (Edge &value: original_edges)
    {
        this->edges[value.get_edge_id()] = &value;
    }

    this->create_adjacency_list();
}


DynamicGraph::DynamicGraph(std::map<int, Vertex *> vertices, std::map<int, Edge *> edges,
                           std::map<std::string, AttributeType> vertex_schema,
                           std::map<std::string, AttributeType> edge_schema,
                           bool undirected, std::map<int, std::vector<Edge *>> adjacency_list) :
    vertices(vertices),
    edges(edges),
    vertex_schema(vertex_schema),
    edge_schema(edge_schema),
    undirected(undirected),
    adjacency_list(adjacency_list)
{
}


void DynamicGraph::create_adjacency_list()
{
    for (auto const &id_vertex_pair : this->vertices) {
        adjacency_list[id_vertex_pair.second->get_vertex_id()] = std::vector<Edge *>();
    }
    for (auto const &id_edge_pair : this->edges) {
        adjacency_list[id_edge_pair.second->get_from_vertex_id()].push_back(id_edge_pair.second);
    }
}


const std::map<int, Vertex *> & DynamicGraph::get_vertices()
{
    return this->vertices;
}
const std::map<int, Edge *> & DynamicGraph::get_edges()
{
    return this->edges;
}
const std::map<std::string, AttributeType> & DynamicGraph::get_vertex_schema()
{
    return this->vertex_schema;
}
const std::map<std::string, AttributeType> & DynamicGraph::get_edge_schema()
{
    return this->edge_schema;
}

bool DynamicGraph::is_undirected()
{
    return this->undirected;
}

const std::map<int, std::vector<Edge *>> & DynamicGraph::get_adjacency_list()
{
    return this->adjacency_list;
}



DynamicGraph DynamicGraph::create_subgraph_with_bounded_edges(timestamp_t upper_time_limit, timestamp_t lower_time_limit)
{
    std::map<int, Edge *> new_edges;

    for (auto const &id_edge_pair : this->edges) {
        if (id_edge_pair.second->get_timestamp() >= lower_time_limit
            && id_edge_pair.second->get_timestamp() < upper_time_limit) {
            new_edges[id_edge_pair.first] = id_edge_pair.second;
        }
    }
    std::map<int, std::vector<Edge *>> new_adjacency_list;

    for (auto &id_list_pair : this->adjacency_list)
    {
        std::vector<Edge *> new_edge_list;
        for (auto &edge : id_list_pair.second)
        {
            if (edge->get_timestamp() >= lower_time_limit
                && edge->get_timestamp() < upper_time_limit) {
                new_edge_list.push_back(edge);
            }
        }
        new_adjacency_list[id_list_pair.first] = new_edge_list;
    }

    return DynamicGraph(vertices, new_edges, vertex_schema, edge_schema, undirected, new_adjacency_list);
}



std::vector<DynamicGraph> DynamicGraph::create_subgraph_instances(std::vector<timestamp_t> positive_event_times,
                                                                  std::vector<timestamp_t> negative_event_times,
                                                                  timestamp_t primary_time_unit)
{
    std::vector<DynamicGraph> dynamic_graphs;
    for (auto positive_event_time : positive_event_times) {
        DynamicGraph dg = this->create_subgraph_with_bounded_edges(positive_event_time,
                                                                     positive_event_time - (5 * primary_time_unit));
        dynamic_graphs.push_back(dg);
    }
    for (auto negative_event_time : negative_event_times) {
        DynamicGraph dg = this->create_subgraph_with_bounded_edges(negative_event_time,
                                                                   negative_event_time - (5 * primary_time_unit));
        dynamic_graphs.push_back(dg);
    }
    return dynamic_graphs;
}


void DynamicGraph::identify_edge_events(std::string attr_name, int label, int sample_size,
                                        RandomGenerator * random_generator, std::vector<int> & edge_ids,
                                        std::vector<std::vector<int>> & vertex_ids,
                                        std::vector<timestamp_t> & vertex_timestamps)
{

    std::vector<int> result_edge_ids;
    std::vector<int> from_vertex_ids;
    std::vector<timestamp_t> edge_timestamps;
    std::vector<int> to_vertex_ids;

    for (auto &kv : edges)
    {
        int edge_attr_val = (int) kv.second->get_attributes().at(attr_name);
        if (edge_attr_val  == label && kv.second->get_edge_id() == kv.second->get_original_edge_id())
        {
            result_edge_ids.push_back(kv.first);
            from_vertex_ids.push_back(kv.second->get_from_vertex_id());
            edge_timestamps.push_back(kv.second->get_timestamp());
            to_vertex_ids.push_back(kv.second->get_to_vertex_id());
        }
    }


    std::vector<int> selected_ids = random_generator->generate_random_int_vector(sample_size, result_edge_ids.size());


    std::vector<int> final_from_vertex_ids;
    std::vector<int> final_to_vertex_ids;
    for (auto &i : selected_ids)
    {
        edge_ids.push_back(result_edge_ids[i]);
        final_from_vertex_ids.push_back(from_vertex_ids[i]);
        final_to_vertex_ids.push_back(to_vertex_ids[i]);
        vertex_timestamps.push_back(edge_timestamps[i]);
    }
    vertex_ids.push_back(final_from_vertex_ids);
    vertex_ids.push_back(final_to_vertex_ids);
}


//
// Created by karbal on 2.1.18.
//

#ifndef EWALDIS_CPP_OCCUPIEDGRAPH_H
#define EWALDIS_CPP_OCCUPIEDGRAPH_H


#include <set>
#include "DynamicGraph.h"

class OccupiedGraph
{
private:
    std::vector<DynamicGraph> & original_graph_instances;
    const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_positive;
    // here we have (instance_id, edge_id) pairs of already selected edges:
    std::set<std::pair<int, int>> all_occupied_original_edges_ids_with_instances;
    // here we have (instance_id, edge_id) pairs of edge candidates:
    std::set<std::pair<int, int>> original_edges_ids_with_instances_in_combinations;
    // how many positive instances we are working with:
    int n_positive;

    // here we keep the occupied vertices and also their vertex_set id
    // mapping (instance_id, vertex_id) -> vertex_set_id
    std::map<int, std::map<int, int>> vertex_mapping;
    // mapping (instance_id, vertex_set_id) -> vertex_id
    std::map<int, std::map<int, int>> vertex_mapping_reversed;
    // counter of vertex sets
    int vertex_sets_current_id;


    // combinations are elements of dictionary, in which key is pair of vertices (FROM, TO)
    // and value is the combination FROM and TO are ids of vertex sets but TO also may be None.
    // each combination is a list in which there is an element for each instance
    // for each such instance, there is a set of edges.
    std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> combinations;
    // here we have the following mapping (instance_id, to_vertex_id) -> (from_vertex_set_id, to_vertex_set_id)
    // and it is when moving to_vertex_id from None to a specific vertex set
    std::map<std::pair<int, int>, std::set<std::pair<int, int>>> to_vertex_id_keys;





    void add_to_combinations(std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & combinations,
                             std::pair<int, int> from_to_id_pair,
                             int instance_id, Edge * edge, int n_positive);


protected:
    virtual void prepare_initial_edge_combination(std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & combinations,
                                                  std::map<std::pair<int, int>, std::set<std::pair<int, int>>> & to_vertex_id_keys);

public:
    ~OccupiedGraph();
    OccupiedGraph(std::vector <DynamicGraph> & graph_instances,
                  const std::vector <std::vector<int>> & positive_event_vertices,
                  const std::vector <std::vector<int>> & positive_event_edges,
                  const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_positive);



    const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & get_edge_combinations();

    virtual int get_vertex_set(int instance_id, int vertex_id);

    void occupy_edges(const std::vector<Edge*> & edges_to_be_occupied);
};


#endif //EWALDIS_CPP_OCCUPIEDGRAPH_H

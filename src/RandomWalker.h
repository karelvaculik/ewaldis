//
// Created by karbal on 29.12.17.
//

#ifndef EWALDIS_CPP_RANDOMWALKER_H
#define EWALDIS_CPP_RANDOMWALKER_H

#include <math.h>
#include <algorithm>
#include <vector>
#include <set>
#include "DynamicGraph.h"
#include "Suitabilities.h"
#include "commonutils.h"
#include "RandomGenerator.h"

class RandomWalker {
private:
    bool use_vertex_attributes;
    timestamp_t time_unit_primary;
    timestamp_t time_unit_secondary;
    int random_walks;
    double prob_restart;
    RandomGenerator * random_generator;

    bool prob_succeed(double prob);

    double expon_pdf(double x, double lamb, timestamp_t time_unit);


    std::pair<void *, double> select_primary_edge(DynamicGraph * graph, int current_node_id, timestamp_t start_time,
                                                  std::set<int> already_visited_edges);

    std::pair<void *, double>
    select_secondary_edge(DynamicGraph * graph, int current_node_id, timestamp_t secondary_event_start_time,
                          timestamp_t secondary_edge_expected_timestamp, Edge * primary_edge,
                          std::set<int> already_visited_edges);

    void perform_one_walk(std::vector<DynamicGraph> & graph_instances, std::vector<int> & current_nodes,
//    void perform_one_walk(DynamicGraph * graph, std::vector<int> & current_nodes,
                          std::vector<timestamp_t > & start_nodes_times, int n_positive_instances,
                          std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_positive,
                          std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_negative,
                          std::map<std::pair<int, int>, double> &edge_priors,
                          std::vector<std::set<int>> & already_visited_edges);

    Suitabilities prepare_scores(std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & pos,
                                 std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & neg,
                                 std::map<std::pair<int, int>, double> & priors);


public:
    RandomWalker(bool use_vertex_attributes,
                 timestamp_t time_unit_primary,
                 timestamp_t time_unit_secondary,
                 int random_walks,
                 double prob_restart,
                 RandomGenerator * random_generator);

    double edge_similarity(DynamicGraph * graph, Edge * primary_edge, Edge * secondary_edge,
                           timestamp_t secondary_event_start_time,
                           timestamp_t secondary_edge_expected_timestamp);

    Suitabilities compute_suitabilities(std::vector<DynamicGraph> & graph_instances,
//    Suitabilities compute_suitabilities(DynamicGraph * graph,
                                        std::vector<std::vector<int>> &positive_event_vertices,
                                        std::vector<timestamp_t> &positive_event_times,
                                        std::vector<std::vector<int>> &positive_event_edges,
                                        std::vector<std::vector<int>> &negative_event_vertices,
                                        std::vector<timestamp_t> &negative_event_times,
                                        std::vector<std::vector<int>> &negative_event_edges);
};


#endif //EWALDIS_CPP_RANDOMWALKER_H

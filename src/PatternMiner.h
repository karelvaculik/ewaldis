//
// Created by karbal on 29.12.17.
//

#ifndef EWALDIS_CPP_PATTERNMINER_H
#define EWALDIS_CPP_PATTERNMINER_H


#include <vector>
#include "DynamicGraph.h"
#include "commonutils.h"
#include "Pattern.h"
#include "RandomGenerator.h"

class PatternMiner
{
private:
    std::vector<DynamicGraph> & graph_instances;
//    DynamicGraph * graph;
    std::vector<std::vector<int>> & positive_event_vertices;
    std::vector<timestamp_t> & positive_event_times;
    std::vector<std::vector<int>> & positive_event_edges;
    std::vector<std::vector<int>> & negative_event_vertices;
    std::vector<timestamp_t> & negative_event_times;
    std::vector<std::vector<int>> & negative_event_edges;
    bool use_vertex_attributes;
    timestamp_t time_unit_primary;
    timestamp_t time_unit_secondary;
    int random_walks;
    double prob_restart;
    int max_pattern_edges;
    int evolution_epochs;
    int evolution_subepochs;
    RandomGenerator * random_generator;
    bool use_simple_init;
    bool use_uniform_crossover;
    bool limit_negative_population;

public:
    PatternMiner(std::vector<DynamicGraph> & graph_instances,
//    PatternMiner(DynamicGraph * graph,
                 std::vector<std::vector<int>> & positive_event_vertices,
                 std::vector<timestamp_t> & positive_event_times,
                 std::vector<std::vector<int>> & positive_event_edges,
                 std::vector<std::vector<int>> & negative_event_vertices,
                 std::vector<timestamp_t> & negative_event_times,
                 std::vector<std::vector<int>> & negative_event_edges,
                 bool use_vertex_attributes,
                 timestamp_t time_unit_primary,
                 timestamp_t time_unit_secondary,
                 int random_walks,
                 double prob_restart,
                 int max_pattern_edges,
                 int evolution_epochs,
                 int evolution_subepochs,
                 RandomGenerator * random_generator,
                 bool use_simple_init,
                 bool use_uniform_crossover,
                 bool limit_negative_population
    );
    Pattern mine_pattern(std::vector<std::vector<double>> &populations_fitness,
                             std::vector<std::vector<double>> &negative_populations_fitness, bool verbose = false);
    std::vector<std::vector<double>> evaluate_pattern(Pattern * pattern,
                                                      std::vector<DynamicGraph> & evaluator_graph_instances,
//                                                      DynamicGraph * graph,
                                                      std::vector<std::vector<int>> & positive_vertex_ids_test,
                                                      std::vector<timestamp_t> & positive_starting_times_test,
                                                      std::vector<std::vector<int>> & negative_vertex_ids_test,
                                                      std::vector<timestamp_t> & negative_starting_times_test,
                                                      int random_walks);

};


#endif //EWALDIS_CPP_PATTERNMINER_H

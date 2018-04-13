//
// Created by karbal on 29.12.17.
//

#include <iostream>
#include "PatternMiner.h"
#include "RandomWalker.h"
#include "GeneticPatternExtractor.h"
#include "PatternChecker.h"
#include "commonutilstemplated.h"


PatternMiner::PatternMiner(std::vector<DynamicGraph> & graph_instances,
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
)
        : graph_instances(graph_instances),
          positive_event_vertices(positive_event_vertices),
          positive_event_times(positive_event_times),
          positive_event_edges(positive_event_edges),
          negative_event_vertices(negative_event_vertices),
          negative_event_times(negative_event_times),
          negative_event_edges(negative_event_edges),
          use_vertex_attributes(use_vertex_attributes),
          time_unit_primary(time_unit_primary),
          time_unit_secondary(time_unit_secondary),
          random_walks(random_walks),
          prob_restart(prob_restart),
          max_pattern_edges(max_pattern_edges),
          evolution_epochs(evolution_epochs),
          evolution_subepochs(evolution_subepochs),
          random_generator(random_generator),
          use_simple_init(use_simple_init),
          use_uniform_crossover(use_uniform_crossover),
          limit_negative_population(limit_negative_population)
{
}


Pattern PatternMiner::mine_pattern(std::vector<std::vector<double>> &populations_fitness,
                                   std::vector<std::vector<double>> &negative_populations_fitness, bool verbose)
{
    if (verbose) println("[STARTED] RANDOM WALKS");

    RandomWalker random_walker = RandomWalker(use_vertex_attributes, time_unit_primary,
                                              time_unit_secondary, random_walks, prob_restart, random_generator);

    Suitabilities suitabilities = random_walker.compute_suitabilities(graph_instances, positive_event_vertices,
                                                                      positive_event_times,
                                                                      positive_event_edges, negative_event_vertices,
                                                                      negative_event_times, negative_event_edges);

    if (verbose) println("[FINISHED] RANDOM WALKS");
    if (verbose) println("[STARTED] PATTERN EXTRACTION");

    GeneticPatternExtractor genetic_pattern_extractor = GeneticPatternExtractor(max_pattern_edges, evolution_epochs,
                                                                                evolution_subepochs, random_generator,
                                                                                use_simple_init, use_uniform_crossover,
                                                                                limit_negative_population);
    Pattern pattern = genetic_pattern_extractor.extract_pattern(graph_instances,
                                                                positive_event_vertices,
                                                                positive_event_times,
                                                                positive_event_edges,
                                                                negative_event_vertices,
                                                                negative_event_edges,
                                                                &suitabilities, populations_fitness,
                                                                negative_populations_fitness,
                                                                verbose);
    if (verbose) println("[FINISHED] PATTERN EXTRACTION");
    return pattern;
}


std::vector<std::vector<double>> PatternMiner::evaluate_pattern(Pattern * pattern,
                                                                std::vector<DynamicGraph> & evaluator_graph_instances,
                                                                std::vector<std::vector<int>> & positive_vertex_ids_test,
                                                                std::vector<timestamp_t> & positive_starting_times_test,
                                                                std::vector<std::vector<int>> & negative_vertex_ids_test,
                                                                std::vector<timestamp_t> & negative_starting_times_test,
                                                                int random_walks)
{
    std::vector<double> positive_scores;
    std::vector<double> negative_scores;


    std::vector<double> edge_weights;
    double sum_of_scores = 0.0;
    for (auto &score : pattern->get_scores())
    {
        edge_weights.push_back(score);
        sum_of_scores += score;
    }
    for (int i = 0; i < edge_weights.size(); ++i)
    {
        edge_weights[i] = edge_weights.at(i) / sum_of_scores;
    }

    const unsigned long n_positive = positive_vertex_ids_test.at(0).size();
    const unsigned long n_negative = negative_vertex_ids_test.at(0).size();

    double score;
    PatternChecker pattern_checker = PatternChecker(use_vertex_attributes,
                                                    time_unit_secondary,
                                                    random_walks,
                                                    random_generator);

    // get the tranposed vectors of vectors
    std::vector<std::vector<int>> transposed_positive_vertex_ids_test(positive_vertex_ids_test[0].size(), vector<int>(positive_vertex_ids_test.size()));
    std::vector<std::vector<int>> transposed_negative_vertex_ids_test(negative_vertex_ids_test[0].size(), vector<int>(negative_vertex_ids_test.size()));
    for (size_t i = 0; i < positive_vertex_ids_test.size(); ++i)
        for (size_t j = 0; j < positive_vertex_ids_test[0].size(); ++j)
            transposed_positive_vertex_ids_test[j][i] = positive_vertex_ids_test[i][j];
    for (size_t i = 0; i < negative_vertex_ids_test.size(); ++i)
        for (size_t j = 0; j < negative_vertex_ids_test[0].size(); ++j)
            transposed_negative_vertex_ids_test[j][i] = negative_vertex_ids_test[i][j];

    for (int j = 0; j < n_positive; ++j)
    {

        score = pattern_checker.check_pattern_in_instance(// &graph_instances.at(j),
                                                          &evaluator_graph_instances.at(j),
                                                          pattern, transposed_positive_vertex_ids_test.at(j),
                                                          positive_starting_times_test.at(j),
                                                          edge_weights);

        positive_scores.push_back(score);
    }

    for (int k = 0; k < n_negative; ++k)
    {
        score = pattern_checker.check_pattern_in_instance(// &graph_instances.at(k),
                                                          &evaluator_graph_instances.at(k),
                                                          pattern, transposed_negative_vertex_ids_test.at(k),
                                                          negative_starting_times_test.at(k),
                                                          edge_weights);
        negative_scores.push_back(score);
    }

    std::vector<std::vector<double>> all_scores = {positive_scores, negative_scores};

    return all_scores;
}


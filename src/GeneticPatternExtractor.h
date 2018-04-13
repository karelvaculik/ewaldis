//
// Created by karel on 2.1.18.
//

#ifndef EWALDIS_CPP_PATTERNEXTRACTOR_H
#define EWALDIS_CPP_PATTERNEXTRACTOR_H


#include <fstream>
#include <math.h>
#include <algorithm>
#include "DynamicGraph.h"
#include "Suitabilities.h"
#include "PatternEdge.h"
#include "OccupiedGraph.h"
#include "RandomGenerator.h"
#include "Pattern.h"

class GeneticPatternExtractor
{
private:
    int max_pattern_edges;
    int evolution_epochs;
    int evolution_subepochs;
    RandomGenerator * random_generator;
    bool use_simple_init;
    bool use_uniform_crossover;

    bool limit_negative_population;


    std::optional<PatternEdge> extract_pattern_edge(OccupiedGraph *occupied_graph, std::vector<DynamicGraph> &graph_instances,
                                                        Suitabilities *suitabilities, int n_positive, int n_negative,
                                                        std::vector<std::vector<double>> &populations_fitness,
                                                        std::vector<std::vector<double>> &negative_populations_fitness,
                                                        int current_pattern_edge_index);

    void initialize_population(const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
                               int n_positive, std::vector<std::vector<std::vector<Edge *>>> &population_positive,
                               std::vector<std::pair<int, int>> &population_keys,
                               std::vector<std::vector<double>> &population_fitness,
                               std::vector<int> &subpopulation_widths);

    void initialize_population_new(const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
                               int n_positive, std::vector<std::vector<std::vector<Edge *>>> &population_positive,
                               std::vector<std::pair<int, int>> &population_keys,
                               std::vector<std::vector<double>> &population_fitness,
                               std::vector<int> &subpopulation_widths);

    void apply_crossover_operator_to_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                                                std::vector<std::vector<double>> &current_fitness);


    void crossover_operator(std::vector<Edge *> &parent_1, std::vector<Edge *> &parent_2,
                            std::vector<Edge *> &child_1, std::vector<Edge *> &child_2);

    void apply_mutation_operator_to_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                                               std::vector<std::pair<int, int>> &population_keys,
                                               const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
                                               std::vector<std::vector<double>> &current_fitness);

    void mutation_operator(std::vector<Edge *> &old_individuum,
                           const std::vector<std::vector<Edge * >> &possible_combinations,
                           std::vector<Edge *> &new_individuum);


    std::vector<double> occupy_best_negative_edges(std::vector<DynamicGraph> &graph_instances,
                                                       std::vector<std::vector<std::vector<Edge *>>> &positive_population,
                                                       int n_positive, int n_negative, Suitabilities *suitabilities,
                                                       int n_subepochs,
                                                       std::vector<std::vector<double>> &negative_populations_fitness,
                                                       int current_pattern_edge_index, int current_epoch);

    void compute_fitness(std::vector<std::vector<std::vector<Edge *>>> &population_positive,
                         std::vector<std::vector<double>> &positive_population_fitness,
                         std::vector<double> &negative_individual_fitness,
                         Suitabilities *suitabilities);

    void select_new_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                               std::vector<std::vector<double>> &fitness,
                               std::vector<int> &subpopulation_widths);


    PatternEdge select_best_individual(std::vector<std::vector<std::vector<Edge *>>> &population,
                                       std::vector<std::vector<double>> &fitness);



    double compute_positive_fitness_of_individual(std::vector<Edge *> &individual,
                                                  Suitabilities *suitabilities);

    std::vector<std::vector<std::vector<Edge *>>> prepare_negative_combinations(
            std::vector<DynamicGraph> &original_graph_instances,
            std::vector<std::vector<std::vector<Edge *>>> &positive_population,
            Suitabilities *suitabilities,
            int n_positive, int n_negative);

    void prepare_negative_population(std::vector<std::vector<std::vector<Edge *>>> &positive_population,
                                     std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
                                     Suitabilities *suitabilities,
                                     int n_positive, int n_negative,
                                     std::vector<std::vector<std::vector<Edge *>>> &negative_population,
                                     std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
                                     std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
                                     std::vector<int> &subpopulation_widths);

    void prepare_negative_population_new(std::vector<std::vector<std::vector<Edge *>>> &positive_population,
                                     std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
                                     Suitabilities *suitabilities,
                                     int n_positive, int n_negative,
                                     std::vector<std::vector<std::vector<Edge *>>> &negative_population,
                                     std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
                                     std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
                                     std::vector<int> &subpopulation_widths);

    void
    compute_partial_negative_fitness_for_individual(std::vector<Edge *> &positive_part,
                                                    std::vector<Edge *> &negative_part,
                                                    Suitabilities *suitabilities,
                                                    int n_positive,
                                                    std::vector<double> &partial_fitness,
                                                    std::vector<int> &partial_counts);

    void apply_crossover_operator_to_negative_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                                                         std::vector<std::vector<std::vector<double>>> &current_partial_fitnesses,
                                                         std::vector<std::vector<std::vector<int>>> &current_partial_counts);

    void apply_mutation_operator_to_negative_population(
            std::vector<std::vector<std::vector<Edge *>>> &negative_population,
            std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
            std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
            std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
            std::vector<std::vector<std::vector<Edge *>>> &positive_population,
            Suitabilities *suitabilities,
            int n_positive);

    std::vector<std::vector<double>>
    compute_negative_fitness(std::vector<std::vector<std::vector<double>>> & current_partial_fitnesses,
                             std::vector<std::vector<std::vector<int>>> & current_partial_counts);

    void select_new_negative_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                                        std::vector<std::vector<double>> &fitness,
                                        std::vector<int> &subpopulation_widths,
                                        std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
                                        std::vector<std::vector<std::vector<int>>> &negative_partial_counts);

    std::vector<double>
    select_best_negative_individual_for_each_positive(
            std::vector<std::vector<std::vector<Edge *>>> &negative_population,
            std::vector<std::vector<double>> &negative_fitness);

    void crossover_operator_negative(std::vector<Edge *> &parent_1, std::vector<Edge *> &parent_2,
                                     std::vector<double> &parent_fit_1, std::vector<double> &parent_fit_2,
                                     std::vector<int> &parent_cnt_1, std::vector<int> &parent_cnt_2,
                                     std::vector<Edge *> &child_1, std::vector<Edge *> &child_2,
                                     std::vector<double> &child_fit_1, std::vector<double> &child_fit_2,
                                     std::vector<int> &child_cnt_1, std::vector<int> &child_cnt_2);

    void
    mutation_operator_negative(std::vector<Edge *> &old_individuum,
                               std::vector<std::vector<Edge *>> &possible_combinations,
                               std::vector<double> &part_fit, std::vector<int> &part_cnt,
                               std::vector<Edge *> &positive_individual,
                               Suitabilities *suitabilities,
                               int n_positive, std::vector<Edge *> &new_individuum,
                               std::vector<double> &new_part_fitness, std::vector<int> &new_part_counts);

    std::pair<int, double> mutated_negative_edge_score(
            Suitabilities *suitabilities,
            int mutation_point, int n_positive, int new_edge_id, std::vector<Edge *> &positive_individual);

public:

    GeneticPatternExtractor(int max_pattern_edges, int evolution_epochs, int evolution_subepochs,
                            RandomGenerator * random_generator,
                            bool use_simple_init,
                            bool use_uniform_crossover,
                            bool limit_negative_population);

    Pattern extract_pattern(std::vector<DynamicGraph> & graph_instances,
                         std::vector<std::vector<int>> & positive_event_vertices,
                         std::vector<timestamp_t> & positive_event_times,
                         std::vector<std::vector<int>> & positive_event_edges,
                         std::vector<std::vector<int>> & negative_event_vertices,
                         std::vector<std::vector<int>> & negative_event_edges,
                         Suitabilities * suitabilities,
                            std::vector<std::vector<double>> & populations_fitness,
                            std::vector<std::vector<double>> & negative_populations_fitness,
                            bool verbose);

};


#endif //EWALDIS_CPP_PATTERNEXTRACTOR_H

//
// Created by karel on 2.1.18.
//

#include <ctime>
#include <fstream>
#include <iostream>
#include <chrono>
#include "GeneticPatternExtractor.h"
#include "Pattern.h"
#include "commonutilstemplated.h"


GeneticPatternExtractor::GeneticPatternExtractor(int max_pattern_edges, int evolution_epochs, int evolution_subepochs,
                                                 RandomGenerator * random_generator,
                                                 bool use_simple_init,
                                                 bool use_uniform_crossover,
                                                 bool limit_negative_population) : max_pattern_edges(max_pattern_edges),
                                                                                     evolution_epochs(evolution_epochs),
                                                                                     evolution_subepochs(evolution_subepochs),
                                                                                     random_generator(random_generator),
                                                                               use_simple_init(use_simple_init),
                                                                               use_uniform_crossover(use_uniform_crossover),
                                                                                   limit_negative_population(limit_negative_population)
{
}


Pattern GeneticPatternExtractor::extract_pattern(std::vector<DynamicGraph> & graph_instances,
//Pattern GeneticPatternExtractor::extract_pattern(DynamicGraph * graph,
                                              std::vector<std::vector<int>> & positive_event_vertices,
                                              std::vector<timestamp_t> & positive_event_times,
                                              std::vector<std::vector<int>> & positive_event_edges,
                                              std::vector<std::vector<int>> & negative_event_vertices,
                                              std::vector<std::vector<int>> & negative_event_edges,
                                              Suitabilities * suitabilities,
                                                 std::vector<std::vector<double>> & populations_fitness,
                                                 std::vector<std::vector<double>> & negative_populations_fitness,
                                                 bool verbose)
{
    OccupiedGraph occupied_graph = OccupiedGraph(graph_instances, positive_event_vertices, positive_event_edges,
//    OccupiedGraph occupied_graph = OccupiedGraph(graph, positive_event_vertices, positive_event_edges,
                                                 suitabilities->get_edge_pairs_dictionary_positive());



    int n_positive = positive_event_vertices.at(0).size();
    int n_negative = negative_event_vertices.at(0).size();

    std::vector<PatternEdge> pattern_edges;


    int current_pattern_edge_index = 0;
    while (current_pattern_edge_index < max_pattern_edges) {

        if (verbose) println("Looking for pattern edge no. ", pattern_edges.size()+1);
        auto pattern_edge_or_null = extract_pattern_edge(&occupied_graph, graph_instances, suitabilities, n_positive,
//        auto pattern_edge_or_null = extract_pattern_edge(&occupied_graph, graph, suitabilities, n_positive,
                                                         n_negative, populations_fitness, negative_populations_fitness,
                                                         current_pattern_edge_index);

        if (!pattern_edge_or_null.has_value())
        {
            break;
        }
        PatternEdge pattern_edge = pattern_edge_or_null.value();
        if (pattern_edge.get_pattern_score() < 0.1)
        {
            break;
        }


        pattern_edges.push_back(pattern_edge);

        occupied_graph.occupy_edges(pattern_edge.get_pattern_edge());
        current_pattern_edge_index++;
    }

    // no create the Pattern from all necessary information:
    std::vector<std::pair<int, int>> pattern_vertex_pairs;
    std::vector<bool> directions;
    std::vector<std::vector<timestamp_t>> timestamps;
    std::vector<std::vector<std::map<std::string, double>>> attributes;
    std::vector<double> scores;

    for (auto& pattern_edge : pattern_edges)
    {
        auto pair = std::make_pair(0, 0);
        bool pair_found = false;
        bool pair_direction = false;
        std::vector<timestamp_t> pair_timestamps;
        std::vector<std::map<std::string, double>> pair_attributes;
        for (int instance_id = 0; instance_id < pattern_edge.get_pattern_edge().size(); ++instance_id) {
            Edge * edge = pattern_edge.get_pattern_edge().at(instance_id);
            if (edge != nullptr)
            {
                if (!pair_found)
                {
                    pair.first = occupied_graph.get_vertex_set(instance_id, edge->get_from_vertex_id());
                    pair.second = occupied_graph.get_vertex_set(instance_id, edge->get_to_vertex_id());
                    pair_direction = edge->get_direction();
                    pair_found = true;
                }
                pair_timestamps.push_back(positive_event_times.at(instance_id) - edge->get_timestamp());
                pair_attributes.push_back(edge->get_attributes());
            }
            // if there is no edge for this instance, push back empty
            else
            {
                pair_timestamps.push_back(0.0);
                pair_attributes.push_back(std::map<std::string, double>());
            }
        }
        pattern_vertex_pairs.push_back(pair);
        directions.push_back(pair_direction);
        timestamps.push_back(pair_timestamps);
        attributes.push_back(pair_attributes);
        scores.push_back(pattern_edge.get_pattern_score());
    }

    Pattern pattern = Pattern(pattern_vertex_pairs, directions,
                              timestamps, attributes, scores);

    return pattern;
}


std::optional<PatternEdge> GeneticPatternExtractor::extract_pattern_edge(OccupiedGraph *occupied_graph, std::vector<DynamicGraph> &graph_instances,
//std::optional<PatternEdge> GeneticPatternExtractor::extract_pattern_edge(OccupiedGraph *occupied_graph, DynamicGraph * graph,
                                                                         Suitabilities *suitabilities, int n_positive, int n_negative,
                                                                         std::vector<std::vector<double>> &populations_fitness,
                                                                         std::vector<std::vector<double>> &negative_populations_fitness,
                                                                         int current_pattern_edge_index)
{
    // computes the whole genetic algorithm in order to find a pattern edge
    std::vector<std::vector<std::vector<Edge *>>> population_positive;
    std::vector<std::pair<int, int>> population_keys;
    std::vector<std::vector<double>> population_fitness;
    std::vector<int> subpopulation_widths;

    if (use_simple_init)
    {
        initialize_population(occupied_graph->get_edge_combinations(), n_positive,
                              population_positive, population_keys, population_fitness, subpopulation_widths);
    }
    else
    {
        initialize_population_new(occupied_graph->get_edge_combinations(), n_positive,
                                  population_positive, population_keys, population_fitness, subpopulation_widths);
    }


    // check whether there are some individuals in the population
    bool empty_population = true;
    for (auto &subpopulation : population_positive) {
        if (subpopulation.size() > 0) {
            empty_population = false;
            break;
        }
    }
    if (empty_population) {
        return {};
    }

    for (int i = 0; i < evolution_epochs; ++i) {
        apply_crossover_operator_to_population(population_positive, population_fitness);
        apply_mutation_operator_to_population(population_positive, population_keys,
                                              occupied_graph->get_edge_combinations(), population_fitness);

        std::vector<double> negative_fitness = occupy_best_negative_edges(graph_instances, population_positive,
//        std::vector<double> negative_fitness = occupy_best_negative_edges(graph, population_positive,
                                                                          n_positive, n_negative, suitabilities,
                                                                          evolution_subepochs,
                                                                          negative_populations_fitness, current_pattern_edge_index, i);

        compute_fitness(population_positive, population_fitness, negative_fitness, suitabilities);

        // store the fitness values
        std::vector<double> population_fitness_flattened;
        for (auto & x : population_fitness)
        {
            for (auto & y : x)
            {
                population_fitness_flattened.push_back(y);
            }
        }
        std::vector<double> quantiles = extract_quantiles(population_fitness_flattened, std::vector<double>{0.0, 0.1, 0.5, 0.9, 1.0});
        populations_fitness[(current_pattern_edge_index * evolution_epochs) + i] = quantiles;

        select_new_population(population_positive, population_fitness, subpopulation_widths);
    }

    return select_best_individual(population_positive, population_fitness);
}


void GeneticPatternExtractor::initialize_population(
        const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
        int n_positive, std::vector<std::vector<std::vector<Edge *>>> &population_positive,
        std::vector<std::pair<int, int>> &population_keys,
        std::vector<std::vector<double>> &population_fitness,
        std::vector<int> &subpopulation_widths)
{
    int xy = 0;
    for (auto &combination_pair : combinations) {
        xy++;
        std::pair<int, int> combination_key = combination_pair.first;
        std::vector<std::vector<Edge *>> combination_list = combination_pair.second;

        // how many edges there are for this key in each instance:
        // (we compute the mean count of the edges across all instances)
        int combination_width = 0;
        for (int i = 0; i < combination_list.size(); ++i) {
            combination_width += combination_list.at(i).size();
        }
        if (combination_width > 0) {
            combination_width = (int) ceil((1.0 * combination_width) / n_positive);
        }

        std::vector<std::vector<Edge *>> subpopulation;
        std::vector<double> subpopulation_fitness;
        if (combination_width > 0) {
            // combination_width tells us how many individuals to prepare
            combination_width = std::max(4, combination_width);


            for (int i = 0; i < combination_width; ++i)
            {
                std::vector<Edge *> individuum;

                for (int j = 0; j < combination_list.size(); ++j) {
                    if (combination_list.at(j).size() == 0)
                    {
                        individuum.push_back(nullptr);
                    }
                    else
                    {
                        int index = random_generator->generate_random_int(0, combination_list.at(j).size());;
                        individuum.push_back(combination_list.at(j).at(index));
                    }
                }

                subpopulation.push_back(individuum);
                subpopulation_fitness.push_back(-1);
            }



        }

        population_positive.push_back(subpopulation);
        population_fitness.push_back(subpopulation_fitness);
        population_keys.push_back(combination_key);
        subpopulation_widths.push_back(combination_width);
    }
}


void GeneticPatternExtractor::initialize_population_new(
        const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
        int n_positive, std::vector<std::vector<std::vector<Edge *>>> &population_positive,
        std::vector<std::pair<int, int>> &population_keys,
        std::vector<std::vector<double>> &population_fitness,
        std::vector<int> &subpopulation_widths)
{
    for (auto &combination_pair : combinations) {
        std::pair<int, int> combination_key = combination_pair.first;
        std::vector<std::vector<Edge *>> combination_list = combination_pair.second;

        // how many edges there are for this key in each instance:
        // (we compute the mean count of the edges across all instances)
        int combination_width = 0;
        for (int i = 0; i < combination_list.size(); ++i) {
            combination_width += combination_list.at(i).size();
        }
        if (combination_width > 0) {
            combination_width = (int) ceil((1.0 * combination_width) / n_positive);
        }

        std::vector<std::vector<Edge *>> subpopulation;
        std::vector<double> subpopulation_fitness;
        if (combination_width > 0) {
            // combination_width tells us how many individuals to prepare
            combination_width = std::max(4, combination_width);

            // here we prepare indices for creating individuals
            // for each instance, there is a vector of indices and each value of this vector, there will be one individual;
            // each inner vector has the length =combination_width or is empty
            std::vector<std::vector<int>> indices_for_individuals;
            for (int k = 0; k < combination_list.size(); ++k)
            {
                // for each instance create a list of indices used for constructing individuals
                if (combination_list.at(k).size() == 0)
                {
                    // for this instance, we cannot use any edges so there are no indices
                    indices_for_individuals.push_back(std::vector<int>());
                }
                else
                {
                    // otherwise find all indices that can be used for building the individuals at this instance
                    // possibly some indices may be repeated
                    indices_for_individuals.push_back(random_generator->generate_random_int_vector(combination_width,
                                                                                                   combination_list.at(k).size()));
                }
            }
            // now use the indices in order to actually build the individuals:
            for (int i = 0; i < combination_width; ++i)
            {
                std::vector<Edge *> individuum;
                for (int j = 0; j < indices_for_individuals.size(); ++j)
                {
                    if (indices_for_individuals.at(j).size() == 0)
                    {
                        individuum.push_back(nullptr);
                    }
                    else
                    {
                        individuum.push_back(combination_list.at(j).at(indices_for_individuals.at(j).at(i)));
                    }
                }
                subpopulation.push_back(individuum);
                subpopulation_fitness.push_back(-1);
            }
        }

        population_positive.push_back(subpopulation);
        population_fitness.push_back(subpopulation_fitness);
        population_keys.push_back(combination_key);
        subpopulation_widths.push_back(combination_width);
    }
}


void GeneticPatternExtractor::apply_crossover_operator_to_population(
        std::vector<std::vector<std::vector<Edge *>>> &population,
        std::vector<std::vector<double>> &current_fitness)
{
    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index) {
        // generate new individuals for each subpopulation that has more than 1 individual
        int subpopulation_size = population.at(subpopulation_index).size();
        if (subpopulation_size > 1) {
            // generate as many new individuals as there are individuals already
            std::vector<std::vector<Edge *>> individuals_to_be_added;
            for (int i = 0; i < (subpopulation_size / 2); ++i) {
                std::pair<int, int> indices = random_generator->generate_random_int_pair(subpopulation_size);
                std::vector<Edge *> new_individuum_1;
                std::vector<Edge *> new_individuum_2;
                crossover_operator(population.at(subpopulation_index).at(indices.first),
                                   population.at(subpopulation_index).at(indices.second),
                                   new_individuum_1, new_individuum_2);
                individuals_to_be_added.push_back(new_individuum_1);
                individuals_to_be_added.push_back(new_individuum_2);
            }
            for (int j = 0; j < individuals_to_be_added.size(); ++j) {
                population.at(subpopulation_index).push_back(individuals_to_be_added.at(j));
                current_fitness.at(subpopulation_index).push_back(-1);
            }
        }
    }
}


void GeneticPatternExtractor::crossover_operator(std::vector<Edge *> &parent_1,
                                                 std::vector<Edge *> &parent_2,
                                                 std::vector<Edge *> &child_1,
                                                 std::vector<Edge *> &child_2)
{
    if (use_uniform_crossover)
    {
        // UNIFORM CROSSOVER
        for (int i = 0; i < parent_1.size(); ++i) {
            if (random_generator->generate_random_double(0.0, 1.0) < 0.5)
            {
                child_1.push_back(parent_1.at(i));
                child_2.push_back(parent_2.at(i));
            }
            else
            {
                child_1.push_back(parent_2.at(i));
                child_2.push_back(parent_1.at(i));
            }
        }
    }
    else
    {
        // SINGLE-POINT CROSSOVER
        int crossing_point = random_generator->generate_random_int(0, parent_1.size() - 1);
        for (int i = 0; i < parent_1.size(); ++i) {
            if (i <= crossing_point) {
                child_1.push_back(parent_1.at(i));
                child_2.push_back(parent_2.at(i));
            } else {
                child_1.push_back(parent_2.at(i));
                child_2.push_back(parent_1.at(i));
            }
        }
    }




}


void GeneticPatternExtractor::apply_mutation_operator_to_population(
        std::vector<std::vector<std::vector<Edge *>>> &population,
        std::vector<std::pair<int, int>> &population_keys,
        const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> &combinations,
        std::vector<std::vector<double>> &current_fitness)
{
    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index) {
        // generate new individuals for each subpopulation that has more than 1 individual
        int subpopulation_size = population.at(subpopulation_index).size();
        if (subpopulation_size > 1) {
            std::vector<std::vector<Edge *>> individuals_to_be_added;
            for (int i = 0; i < subpopulation_size; ++i) {
                int index_1 = random_generator->generate_random_int(0, subpopulation_size);
                std::vector<Edge *> new_individuum;
                mutation_operator(population.at(subpopulation_index).at(index_1), combinations.at(population_keys.at(subpopulation_index)),
                                  new_individuum);

                if (new_individuum.size() > 0) {
                    // if it has size > 0 (i.e. it was really created)
                    individuals_to_be_added.push_back(new_individuum);
                }
            }
            for (int j = 0; j < individuals_to_be_added.size(); ++j) {
                population.at(subpopulation_index).push_back(individuals_to_be_added.at(j));
                current_fitness.at(subpopulation_index).push_back(-1);
            }
        }
    }
}


void GeneticPatternExtractor::mutation_operator(std::vector<Edge *> &old_individuum,
                                                const std::vector<std::vector<Edge *>> &possible_combinations,
                                                std::vector<Edge *> &new_individuum)
{
    int mutation_point = random_generator->generate_random_int(0, old_individuum.size());
    if (old_individuum.at(mutation_point) == nullptr) {
        return;
    }

    // find potential edges that can be used at the mutation point
    int current_edge_original_id = old_individuum.at(mutation_point)->get_original_edge_id();
    std::vector<Edge *> potential_edges;
    for (int i = 0; i < possible_combinations.at(mutation_point).size(); ++i) {
        if (possible_combinations.at(mutation_point).at(i)->get_original_edge_id() != current_edge_original_id) {
            potential_edges.push_back(possible_combinations.at(mutation_point).at(i));
        }
    }
    // if there are some potential edges, pick one (otherwise, don't do anything)
    if (potential_edges.size() > 0)
    {
        for (int j = 0; j < old_individuum.size(); ++j)
        {
            if (j != mutation_point)
            {
                new_individuum.push_back(old_individuum.at(j));
            }
            else
            {
                int new_edge_ind = random_generator->generate_random_int(0, potential_edges.size());
                new_individuum.push_back(potential_edges.at(new_edge_ind));
            }
        }
    }

}


std::vector<double> GeneticPatternExtractor::occupy_best_negative_edges(std::vector<DynamicGraph> &graph_instances,
//std::vector<double> GeneticPatternExtractor::occupy_best_negative_edges(DynamicGraph * graph,
                                                                        std::vector<std::vector<std::vector<Edge *>>> &positive_population,
                                                                        int n_positive, int n_negative, Suitabilities *suitabilities,
                                                                        int n_subepochs,
                                                                        std::vector<std::vector<double>> &negative_populations_fitness,
                                                                        int current_pattern_edge_index, int current_epoch)
{
    // for each positive individual, here are possible combinations:
    std::vector<std::vector<std::vector<Edge *>>> negative_combinations = prepare_negative_combinations(
            graph_instances,
//            graph,
            positive_population,
            suitabilities,
            n_positive,
            n_negative);



    std::vector<std::vector<std::vector<Edge *>>> negative_population;
    std::vector<std::vector<std::vector<double>>> negative_partial_fitnesses;
    std::vector<std::vector<std::vector<int>>> negative_partial_counts;
    std::vector<int> subpopulation_widths;

    if (use_simple_init)
    {
        prepare_negative_population(positive_population, negative_combinations,
                                    suitabilities,
                                    n_positive, n_negative,
                                    negative_population, negative_partial_fitnesses, negative_partial_counts,
                                    subpopulation_widths);
    }
    else
    {
        prepare_negative_population_new(positive_population, negative_combinations,
                                        suitabilities,
                                        n_positive, n_negative,
                                        negative_population, negative_partial_fitnesses, negative_partial_counts,
                                        subpopulation_widths);
    }

    bool empty_population = true;
    for (auto &subpopulation : negative_population) {
        if (subpopulation.size() > 0) {
            empty_population = false;
            break;
        }
    }
    if (empty_population) {
        return std::vector<double>(negative_population.size(), 0.0);
    }

    std::vector<std::vector<double>> negative_fitness;
    for (int i = 0; i < n_subepochs; ++i)
    {
        apply_crossover_operator_to_negative_population(negative_population, negative_partial_fitnesses,
                                                        negative_partial_counts);
        apply_mutation_operator_to_negative_population(negative_population, negative_combinations,
                                                       negative_partial_fitnesses, negative_partial_counts,
                                                       positive_population, suitabilities,
                                                       n_positive);

        negative_fitness = compute_negative_fitness(negative_partial_fitnesses, negative_partial_counts);

        // store the negative fitness values
        std::vector<double> population_fitness_flattened;
        for (auto & x : negative_fitness)
        {
            for (auto & y : x)
            {
                population_fitness_flattened.push_back(y);
            }
        }
        std::vector<double> quantiles = extract_quantiles(population_fitness_flattened, std::vector<double>{0.0, 0.1, 0.5, 0.9, 1.0});
        negative_populations_fitness[(current_pattern_edge_index * evolution_epochs * n_subepochs) + (current_epoch * n_subepochs) + i] = quantiles;

        select_new_negative_population(negative_population, negative_fitness, subpopulation_widths,
                                       negative_partial_fitnesses,
                                       negative_partial_counts);

    }

    std::vector<double> negative_individuals_score = select_best_negative_individual_for_each_positive(
            negative_population,
            negative_fitness);
    return negative_individuals_score;
}

void GeneticPatternExtractor::compute_fitness(std::vector<std::vector<std::vector<Edge *>>> &population_positive,
                                              std::vector<std::vector<double>> &positive_population_fitness,
                                              std::vector<double> &negative_individual_fitness,
//                                                 std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_positive
                                              Suitabilities *suitabilities
)
{
    // clear the vector (we don't cache positive fitness between evolutions at this moment)
    positive_population_fitness.clear();
    int negative_fitness_index = 0;
    for (int subpopulation_index = 0; subpopulation_index < population_positive.size(); ++subpopulation_index)
    {
        std::vector<double> subpopulation_fitness;
        for (int individual_index = 0; individual_index < population_positive.at(subpopulation_index).size(); ++individual_index)
        {
            double positive_score = compute_positive_fitness_of_individual(
                    population_positive.at(subpopulation_index).at(individual_index),
                    suitabilities);
            subpopulation_fitness.push_back(positive_score - negative_individual_fitness.at(negative_fitness_index));
            ++negative_fitness_index;
        }
        positive_population_fitness.push_back(subpopulation_fitness);
    }
}


double GeneticPatternExtractor::compute_positive_fitness_of_individual(std::vector<Edge *> &individual,
//                                                                          std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_positive
                                                                       Suitabilities *suitabilities
)
{
    double score = 0.0;
    int cnt = 0;
    for (int pos_index_1 = 0; pos_index_1 < individual.size() - 1; ++pos_index_1)
    {
        if (individual.at(pos_index_1) == nullptr)
        {
            continue;
        }
        for (int pos_index_2 = pos_index_1 + 1; pos_index_2 < individual.size(); ++pos_index_2)
        {
            if (individual.at(pos_index_2) == nullptr)
            {
                continue;
            }
            std::pair<int, int> pos_key_1 = std::make_pair(pos_index_1, individual.at(pos_index_1)->get_edge_id());
            std::pair<int, int> pos_key_2 = std::make_pair(pos_index_2, individual.at(pos_index_2)->get_edge_id());
            if (suitabilities->get_edge_pairs_dictionary_positive().find(pos_key_1) != suitabilities->get_edge_pairs_dictionary_positive().end())
            {
                if (suitabilities->get_edge_pairs_dictionary_positive().at(pos_key_1).find(pos_key_2) != suitabilities->get_edge_pairs_dictionary_positive().at(pos_key_1).end())
                {
                    score += suitabilities->get_edge_pairs_dictionary_positive().at(pos_key_1).at(pos_key_2);
                    cnt += 1;
                }
            }
        }
    }
    if (cnt > 0)
    {
        return score / cnt;
    }
    else
    {
        return 0.0;
    }
}


void GeneticPatternExtractor::select_new_population(std::vector<std::vector<std::vector<Edge *>>> &population,
                                                    std::vector<std::vector<double>> &fitness,
                                                    std::vector<int> &subpopulation_widths)
{

    std::vector<std::vector<std::vector<Edge *>>> new_population;
    std::vector<std::vector<double>> new_fitness;

    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index)
    {
        std::vector<std::vector<Edge *>> new_subpopulation;
        std::vector<double> new_subpopulation_fitness;
        if (population.at(subpopulation_index).size() > 0)
        {
            double min_fitness = fitness.at(subpopulation_index).at(0);
            for (int i = 0; i < fitness.at(subpopulation_index).size(); ++i) {
                if (fitness.at(subpopulation_index).at(i) < min_fitness)
                {
                    min_fitness = fitness.at(subpopulation_index).at(i);
                }
            }
            std::vector<int> selected_individuals_indices;
            if (min_fitness <= 0.0)
            {
                // if there are some negative fitness values, shift all values so that minimum is zero + small smoothing
                std::vector<double> weights;
                for (auto& fit : fitness.at(subpopulation_index))
                {
                    weights.push_back(fit - min_fitness + 0.001);
                }
                selected_individuals_indices = random_generator->stochastic_universal_sampling(
                        subpopulation_widths.at(subpopulation_index), weights);
            }
            else
            {
                // otherwise just use the fitness values as they are
                selected_individuals_indices = random_generator->stochastic_universal_sampling(
                        subpopulation_widths.at(subpopulation_index), fitness.at(subpopulation_index));
            }

            for (auto& ind : selected_individuals_indices)
            {
                new_subpopulation.push_back(population.at(subpopulation_index).at(ind));
                new_subpopulation_fitness.push_back(fitness.at(subpopulation_index).at(ind));
            }

        }
        new_population.push_back(new_subpopulation);
        new_fitness.push_back(new_subpopulation_fitness);
    }

    population = new_population;
    fitness = new_fitness;

}


PatternEdge GeneticPatternExtractor::select_best_individual(
        std::vector<std::vector<std::vector<Edge *>>> &population,
        std::vector<std::vector<double>> &fitness)
{
    std::vector<Edge*> current_best_individual;
    double current_best_fitness = -1;
    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index)
    {
        for (int individual_index = 0; individual_index < population.at(subpopulation_index).size(); ++individual_index)
        {
            if (fitness.at(subpopulation_index).at(individual_index) > current_best_fitness)
            {
                current_best_fitness = fitness.at(subpopulation_index).at(individual_index);
                current_best_individual = population.at(subpopulation_index).at(individual_index);
            }
        }
    }
    return PatternEdge(current_best_individual, current_best_fitness);
}

std::vector<std::vector<std::vector<Edge *>>> GeneticPatternExtractor::prepare_negative_combinations(
        std::vector<DynamicGraph> & original_graph_instances,
//        DynamicGraph * graph,
        std::vector<std::vector<std::vector<Edge *>>> &positive_population,
        Suitabilities *suitabilities,
        int n_positive, int n_negative)
{
    std::vector<std::vector<std::vector<Edge *>>> negative_combinations;

    for (int subpopulation_index = 0; subpopulation_index < positive_population.size(); ++subpopulation_index)
    {
        for (int individual_index = 0; individual_index < positive_population.at(subpopulation_index).size(); ++individual_index)
        {
            std::vector<std::vector<Edge *>> subcombination;
            for (int i = 0; i < n_negative; ++i)
            {
                subcombination.push_back(std::vector<Edge *>());
            }
            for (int instance_id = 0; instance_id < positive_population.at(subpopulation_index).at(individual_index).size(); ++instance_id)
            {
                if (positive_population.at(subpopulation_index).at(individual_index).at(instance_id) == nullptr)
                {
                    continue;
                }
                std::pair<int, int> id_pair = std::make_pair(instance_id, positive_population.at(subpopulation_index).at(individual_index).at(instance_id)->get_edge_id());
                if (suitabilities->get_edge_pairs_dictionary_negative().find(id_pair) != suitabilities->get_edge_pairs_dictionary_negative().end())
                {
                    std::map<std::pair<int, int>, double> tmp = suitabilities->get_edge_pairs_dictionary_negative().at(id_pair);
                    for (auto & kv : tmp)
                    {
                        int neg_instance_id = kv.first.first;
                        int neg_edge_id = kv.first.second;
                        subcombination.at(neg_instance_id - n_positive).push_back(original_graph_instances.at(neg_instance_id).get_edges().at(neg_edge_id));
//                        subcombination.at(neg_instance_id - n_positive).push_back(graph->get_edges().at(neg_edge_id));
                    }
                }
            }
            negative_combinations.push_back(subcombination);
        }
    }

    return negative_combinations;
}

void
GeneticPatternExtractor::prepare_negative_population(
        std::vector<std::vector<std::vector<Edge *>>> &positive_population,
        std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
        Suitabilities *suitabilities,
        int n_positive, int n_negative,
        std::vector<std::vector<std::vector<Edge *>>> &negative_population,
        std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
        std::vector<int> &subpopulation_widths)
{
    int subcombination_index = 0;
    for (int pos_subpopulation_index = 0; pos_subpopulation_index < positive_population.size(); ++pos_subpopulation_index)
    {
        for (int pos_individual_index = 0; pos_individual_index < positive_population.at(pos_subpopulation_index).size(); ++pos_individual_index)
        {
            // how many edges there are for this key in each instance:
            // (we compute the mean count of the edges across all instances)
            int combination_width = 0;
            for (int i = 0; i < negative_combinations.at(subcombination_index).size(); ++i)
            {
                combination_width += negative_combinations.at(subcombination_index).at(i).size();
            }
            if (combination_width > 0)
            {
                combination_width = (int) ceil((1.0 * combination_width) / negative_combinations.at(subcombination_index).size());
            }
            std::vector<std::vector<Edge *>> subpopulation;
            std::vector<std::vector<double>> subpopulation_fitnesses;
            std::vector<std::vector<int>> subpopulation_counts;
            if (combination_width > 0)
            {
                //  but take minimally 4 and maximally 10:
                if (limit_negative_population)
                {
                    combination_width = std::max(4, std::min(combination_width, 10));
                }
                else
                {
                    combination_width = std::max(4, combination_width);
                }

                for (int i = 0; i < combination_width; ++i)
                {
                    // create a new negative individual:
                    std::vector<Edge *> negative_individual;
                    for (int j = 0; j < negative_combinations.at(subcombination_index).size(); ++j)
                    {
                        if (negative_combinations.at(subcombination_index).at(j).size() == 0)
                        {
                            negative_individual.push_back(nullptr);
                        }
                        else
                        {
                            int index = random_generator->generate_random_int(0, negative_combinations.at(subcombination_index).at(j).size());;
                            negative_individual.push_back(negative_combinations.at(subcombination_index).at(j).at(index));
                        }
                    }
                    // initialise the partial fitnesses and counts to zeros:
                    std::vector<double> partial_fitness(n_negative, 0.0);
                    std::vector<int> partial_counts(n_negative, 0);

                    compute_partial_negative_fitness_for_individual(
                            positive_population.at(pos_subpopulation_index).at(pos_individual_index),
                            negative_individual, suitabilities, n_positive, partial_fitness,
                            partial_counts);

                    subpopulation.push_back(negative_individual);
                    subpopulation_fitnesses.push_back(partial_fitness);
                    subpopulation_counts.push_back(partial_counts);
                }



            }

            negative_population.push_back(subpopulation);
            negative_partial_fitnesses.push_back(subpopulation_fitnesses);
            negative_partial_counts.push_back(subpopulation_counts);
            subpopulation_widths.push_back(combination_width);
            ++subcombination_index;
        }
    }
}


void
GeneticPatternExtractor::prepare_negative_population_new(
        std::vector<std::vector<std::vector<Edge *>>> &positive_population,
        std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
        Suitabilities *suitabilities,
        int n_positive, int n_negative,
        std::vector<std::vector<std::vector<Edge *>>> &negative_population,
        std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
        std::vector<int> &subpopulation_widths)
{
    int subcombination_index = 0;
    for (int pos_subpopulation_index = 0; pos_subpopulation_index < positive_population.size(); ++pos_subpopulation_index)
    {
        for (int pos_individual_index = 0; pos_individual_index < positive_population.at(pos_subpopulation_index).size(); ++pos_individual_index)
        {
            // how many edges there are for this key in each instance:
            // (we compute the mean count of the edges across all instances)
            int combination_width = 0;
            for (int i = 0; i < negative_combinations.at(subcombination_index).size(); ++i)
            {
                combination_width += negative_combinations.at(subcombination_index).at(i).size();
            }
            if (combination_width > 0)
            {
                combination_width = (int) ceil((1.0 * combination_width) / negative_combinations.at(subcombination_index).size());
            }
            std::vector<std::vector<Edge *>> subpopulation;
            std::vector<std::vector<double>> subpopulation_fitnesses;
            std::vector<std::vector<int>> subpopulation_counts;
            if (combination_width > 0)
            {
                //  but take minimally 4 and maximally 10:
                if (limit_negative_population)
                {
                    combination_width = std::max(4, std::min(combination_width, 10));
                }
                else
                {
                    combination_width = std::max(4, combination_width);
                }

                // here we prepare indices for creating individuals
                // for each instance, there is a vector of indices and each value of this vector, there will be one individual;
                // each inner vector has the length =combination_width or is empty
                std::vector<std::vector<int>> indices_for_individuals;
                for (int k = 0; k < negative_combinations.at(subcombination_index).size(); ++k)
                {
                    // for each instance create a list of indices used for constructing individuals
                    if (negative_combinations.at(subcombination_index).at(k).size() == 0)
                    {
                        // for this instance, we cannot use any edges so there are no indices
                        indices_for_individuals.push_back(std::vector<int>());
                    }
                    else
                    {
                        // otherwise find all indices that can be used for building the individuals at this instance
                        // possibly some indices may be repeated
                        indices_for_individuals.push_back(random_generator->generate_random_int_vector(combination_width,
                                                                                                       negative_combinations.at(subcombination_index).at(k).size()));
                    }
                }
                // now use the indices in order to actually build the individuals:
                for (int i = 0; i < combination_width; ++i)
                {
                    std::vector<Edge *> negative_individual;
                    for (int j = 0; j < indices_for_individuals.size(); ++j)
                    {
                        if (indices_for_individuals.at(j).size() == 0)
                        {
                            negative_individual.push_back(nullptr);
                        }
                        else
                        {
                            negative_individual.push_back(negative_combinations.at(subcombination_index).at(j).at(indices_for_individuals.at(j).at(i)));
                        }
                    }
                    // initialise the partial fitnesses and counts to zeros:
                    std::vector<double> partial_fitness(n_negative, 0.0);
                    std::vector<int> partial_counts(n_negative, 0);

                    compute_partial_negative_fitness_for_individual(
                            positive_population.at(pos_subpopulation_index).at(pos_individual_index),
                            negative_individual, suitabilities, n_positive, partial_fitness,
                            partial_counts);

                    subpopulation.push_back(negative_individual);
                    subpopulation_fitnesses.push_back(partial_fitness);
                    subpopulation_counts.push_back(partial_counts);
                }

            }

            negative_population.push_back(subpopulation);
            negative_partial_fitnesses.push_back(subpopulation_fitnesses);
            negative_partial_counts.push_back(subpopulation_counts);
            subpopulation_widths.push_back(combination_width);
            ++subcombination_index;
        }
    }
}




void GeneticPatternExtractor::compute_partial_negative_fitness_for_individual(std::vector<Edge *> &positive_part,
                                                                              std::vector<Edge *> &negative_part,
//                                                                                 std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_negative,
                                                                              Suitabilities *suitabilities,
                                                                              int n_positive,
                                                                              std::vector<double> &partial_fitness,
                                                                              std::vector<int> &partial_counts)
{
    for (int pos_instance_id = 0; pos_instance_id < positive_part.size(); ++pos_instance_id)
    {
        if (positive_part.at(pos_instance_id) == nullptr)
        {
            continue;
        }
        for (int neg_instance_id = 0; neg_instance_id < negative_part.size(); ++neg_instance_id)
        {
            if (negative_part.at(neg_instance_id) == nullptr)
            {
                // if there is None instead of an edge in the individual, skip it
                continue;
            }
            std::pair<int, int> pos_id_pair = std::make_pair(pos_instance_id, positive_part.at(pos_instance_id)->get_edge_id());
            std::pair<int, int> neg_id_pair = std::make_pair(n_positive + neg_instance_id, negative_part.at(neg_instance_id)->get_edge_id());

            if ((suitabilities->get_edge_pairs_dictionary_negative().find(pos_id_pair) != suitabilities->get_edge_pairs_dictionary_negative().end())
                    && (suitabilities->get_edge_pairs_dictionary_negative().at(pos_id_pair).find(neg_id_pair) != suitabilities->get_edge_pairs_dictionary_negative().at(pos_id_pair).end()))
            {
                partial_fitness[neg_instance_id] += suitabilities->get_edge_pairs_dictionary_negative().at(pos_id_pair).at(neg_id_pair);
                partial_counts[neg_instance_id] += 1;
            }
        }
    }
}

void GeneticPatternExtractor::apply_crossover_operator_to_negative_population(
        std::vector<std::vector<std::vector<Edge *>>> &population,
        std::vector<std::vector<std::vector<double>>> &current_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> &current_partial_counts)
{
    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index) {
        // generate new individuals for each subpopulation that has more than 1 individual
        int subpopulation_size = population.at(subpopulation_index).size();
        if (subpopulation_size > 1) {
            // generate as many new individuals as there are individuals already
            std::vector<std::vector<Edge *>> individuals_to_be_added;
            std::vector<std::vector<double>> partial_fitness_values_to_be_added;
            std::vector<std::vector<int>> partial_counts_values_to_be_added;

            for (int i = 0; i < (subpopulation_size / 2); ++i) {
                std::pair<int, int> indices = random_generator->generate_random_int_pair(subpopulation_size);
                std::vector<Edge *> new_individuum_1;
                std::vector<Edge *> new_individuum_2;
                std::vector<double> new_fit_1;
                std::vector<double> new_fit_2;
                std::vector<int> new_cnt_1;
                std::vector<int> new_cnt_2;
                crossover_operator_negative(population.at(subpopulation_index).at(indices.first),
                                            population.at(subpopulation_index).at(indices.second),
                                            current_partial_fitnesses.at(subpopulation_index).at(indices.first),
                                            current_partial_fitnesses.at(subpopulation_index).at(indices.second),
                                            current_partial_counts.at(subpopulation_index).at(indices.first),
                                            current_partial_counts.at(subpopulation_index).at(indices.second),
                                            new_individuum_1, new_individuum_2,
                                            new_fit_1, new_fit_2, new_cnt_1, new_cnt_2);
                individuals_to_be_added.push_back(new_individuum_1);
                individuals_to_be_added.push_back(new_individuum_2);
                partial_fitness_values_to_be_added.push_back(new_fit_1);
                partial_fitness_values_to_be_added.push_back(new_fit_2);
                partial_counts_values_to_be_added.push_back(new_cnt_1);
                partial_counts_values_to_be_added.push_back(new_cnt_2);
            }
            for (int j = 0; j < individuals_to_be_added.size(); ++j) {
                population.at(subpopulation_index).push_back(individuals_to_be_added.at(j));
                current_partial_fitnesses.at(subpopulation_index).push_back(partial_fitness_values_to_be_added.at(j));
                current_partial_counts.at(subpopulation_index).push_back(partial_counts_values_to_be_added.at(j));
            }
        }
    }
}


void GeneticPatternExtractor::crossover_operator_negative(std::vector<Edge *> &parent_1,
                                                          std::vector<Edge *> &parent_2,
                                                          std::vector<double> &parent_fit_1,
                                                          std::vector<double> &parent_fit_2,
                                                          std::vector<int> &parent_cnt_1, std::vector<int> &parent_cnt_2,
                                                          std::vector<Edge *> &child_1,
                                                          std::vector<Edge *> &child_2,
                                                          std::vector<double> &child_fit_1,
                                                          std::vector<double> &child_fit_2,
                                                          std::vector<int> &child_cnt_1, std::vector<int> &child_cnt_2)
{
    if (use_uniform_crossover)
    {

        // UNIFORM CROSSOVER
        for (int i = 0; i < parent_1.size(); ++i) {
            if (random_generator->generate_random_double(0.0, 1.0) < 0.5)
            {
                child_1.push_back(parent_1.at(i));
                child_2.push_back(parent_2.at(i));

                child_fit_1.push_back(parent_fit_1.at(i));
                child_fit_2.push_back(parent_fit_2.at(i));

                child_cnt_1.push_back(parent_cnt_1.at(i));
                child_cnt_2.push_back(parent_cnt_2.at(i));
            }
            else
            {
                child_1.push_back(parent_2.at(i));
                child_2.push_back(parent_1.at(i));

                child_fit_1.push_back(parent_fit_2.at(i));
                child_fit_2.push_back(parent_fit_1.at(i));

                child_cnt_1.push_back(parent_cnt_2.at(i));
                child_cnt_2.push_back(parent_cnt_1.at(i));
            }
        }
    }
    else
    {
        // SINGLE-POINT CROSSOVER
        int crossing_point = random_generator->generate_random_int(0, parent_1.size() - 1);
        for (int i = 0; i < parent_1.size(); ++i) {
            if (i <= crossing_point) {
                child_1.push_back(parent_1.at(i));
                child_2.push_back(parent_2.at(i));

                child_fit_1.push_back(parent_fit_1.at(i));
                child_fit_2.push_back(parent_fit_2.at(i));

                child_cnt_1.push_back(parent_cnt_1.at(i));
                child_cnt_2.push_back(parent_cnt_2.at(i));

            } else {
                child_1.push_back(parent_2.at(i));
                child_2.push_back(parent_1.at(i));

                child_fit_1.push_back(parent_fit_2.at(i));
                child_fit_2.push_back(parent_fit_1.at(i));

                child_cnt_1.push_back(parent_cnt_2.at(i));
                child_cnt_2.push_back(parent_cnt_1.at(i));
            }
        }

    }




}


void GeneticPatternExtractor::apply_mutation_operator_to_negative_population(
        std::vector<std::vector<std::vector<Edge *>>> &negative_population,
        std::vector<std::vector<std::vector<Edge *>>> &negative_combinations,
        std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> &negative_partial_counts,
        std::vector<std::vector<std::vector<Edge *>>> &positive_population,
        Suitabilities *suitabilities,
        int n_positive)
{
    int negative_subpopulation_index = 0;
    for (int pos_subpopulation_index = 0; pos_subpopulation_index < positive_population.size(); ++pos_subpopulation_index)
    {
        for (int pos_individual_index = 0; pos_individual_index < positive_population.at(pos_subpopulation_index).size(); ++pos_individual_index)
        {
            int subpopulation_size = negative_population.at(negative_subpopulation_index).size();
            if (subpopulation_size > 1)
            {
                std::vector<std::vector<Edge *>> new_individuals;
                std::vector<std::vector<double>> new_partial_fitnesses;
                std::vector<std::vector<int>> new_partial_counts;
                for (int i = 0; i < subpopulation_size; ++i)
                {
                    int index_1 = random_generator->generate_random_int(0, subpopulation_size);
                    std::vector<Edge *> new_individuum;
                    std::vector<double> new_partial_fitness;
                    std::vector<int> new_partial_count;

                    mutation_operator_negative(negative_population.at(negative_subpopulation_index).at(index_1),
                                               negative_combinations.at(negative_subpopulation_index),
                                               negative_partial_fitnesses.at(negative_subpopulation_index).at(index_1),
                                               negative_partial_counts.at(negative_subpopulation_index).at(index_1),
                                               positive_population.at(pos_subpopulation_index).at(pos_individual_index),
                                               suitabilities, n_positive,
                                               new_individuum, new_partial_fitness, new_partial_count);

                    if (new_individuum.size() > 0) {
                        // if it has size > 0 (i.e. it was really created)
                        new_individuals.push_back(new_individuum);
                        new_partial_fitnesses.push_back(new_partial_fitness);
                        new_partial_counts.push_back(new_partial_count);
                    }
                }
                for (int j = 0; j < new_individuals.size(); ++j) {
                    negative_population.at(negative_subpopulation_index).push_back(new_individuals.at(j));
                    negative_partial_fitnesses.at(negative_subpopulation_index).push_back(new_partial_fitnesses.at(j));
                    negative_partial_counts.at(negative_subpopulation_index).push_back(new_partial_counts.at(j));

                }

            }

            ++negative_subpopulation_index;
        }
    }
}


void GeneticPatternExtractor::mutation_operator_negative(std::vector<Edge *> &old_individuum,
                                                         std::vector<std::vector<Edge *>> &possible_combinations,
                                                         std::vector<double> &part_fit, std::vector<int> &part_cnt,
                                                         std::vector<Edge *> &positive_individual,
                                                         Suitabilities *suitabilities,
//                                                            std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> &edge_pairs_dictionary_negative,
                                                         int n_positive, std::vector<Edge *> &new_individuum,
                                                         std::vector<double> &new_part_fitness,
                                                         std::vector<int> &new_part_counts)
{
    int mutation_point = random_generator->generate_random_int(0, old_individuum.size());
    if (old_individuum.at(mutation_point) == nullptr) {
        return;
    }

    // find potential edges that can be used at the mutation point
    int current_edge_original_id = old_individuum.at(mutation_point)->get_original_edge_id();
    std::vector<Edge *> potential_edges;
    for (int i = 0; i < possible_combinations.at(mutation_point).size(); ++i) {
        if (possible_combinations.at(mutation_point).at(i)->get_original_edge_id() != current_edge_original_id) {
            potential_edges.push_back(possible_combinations.at(mutation_point).at(i));
        }
    }
    // if there are some potential edges, pick one (otherwise, don't do anything)
    if (potential_edges.size() > 0)
    {
        for (int j = 0; j < old_individuum.size(); ++j)
        {
            if (j != mutation_point)
            {
                new_individuum.push_back(old_individuum.at(j));
                new_part_fitness.push_back(part_fit.at(j));
                new_part_counts.push_back(part_cnt.at(j));
            }
            else
            {
                int new_edge_ind = random_generator->generate_random_int(0, potential_edges.size());
                new_individuum.push_back(potential_edges.at(new_edge_ind));

                std::pair<int, double> cnt_fit_pair = mutated_negative_edge_score(suitabilities,
                                                                                  mutation_point, n_positive,
                                                                                  potential_edges.at(
                                                                                          new_edge_ind)->get_edge_id(),
                                                                                  positive_individual);

                new_part_fitness.push_back(cnt_fit_pair.second);
                new_part_counts.push_back(cnt_fit_pair.first);
            }
        }
    }

}


std::pair<int, double> GeneticPatternExtractor::mutated_negative_edge_score(
        Suitabilities *suitabilities,
        int mutation_point, int n_positive, int new_edge_id, std::vector<Edge *> &positive_individual)
{
    double new_fit = 0.0;
    int new_cnt = 0;
    std::pair<int, int> neg_key = std::make_pair(n_positive + mutation_point, new_edge_id);
    for (int pos_instance_id = 0; pos_instance_id < positive_individual.size(); ++pos_instance_id)
    {
        if (positive_individual.at(pos_instance_id) == nullptr)
        {
            continue;
        }
        std::pair<int, int> pos_key = std::make_pair(pos_instance_id, positive_individual.at(pos_instance_id)->get_edge_id());
        if ((suitabilities->get_edge_pairs_dictionary_negative().find(pos_key) != suitabilities->get_edge_pairs_dictionary_negative().end())
            && (suitabilities->get_edge_pairs_dictionary_negative().at(pos_key).find(neg_key) != suitabilities->get_edge_pairs_dictionary_negative().at(pos_key).end()))
        {
            new_fit += suitabilities->get_edge_pairs_dictionary_negative().at(pos_key).at(neg_key);
            new_cnt += 1;
        }
    }

    return std::make_pair(new_cnt, new_fit);
}



std::vector<std::vector<double>> GeneticPatternExtractor::compute_negative_fitness(
        std::vector<std::vector<std::vector<double>>> & current_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> & current_partial_counts)
{
    std::vector<std::vector<double>> negative_fitness;

    for (int subpopulation_index = 0; subpopulation_index < current_partial_fitnesses.size(); ++subpopulation_index)
    {
        std::vector<double> negative_subpopulation_fitness;
        for (int individual_index = 0; individual_index < current_partial_fitnesses.at(subpopulation_index).size(); ++individual_index)
        {
            double fit_sum = 0.0;
            for (auto& n : current_partial_fitnesses.at(subpopulation_index).at(individual_index))
            {
                fit_sum += n;
            }
            int cnt_sum = 0;
            for (auto& n : current_partial_counts.at(subpopulation_index).at(individual_index))
            {
                cnt_sum += n;
            }
            if (cnt_sum > 0)
            {
                negative_subpopulation_fitness.push_back(fit_sum / cnt_sum);
            }
            else
            {
                negative_subpopulation_fitness.push_back(0.0);
            }
        }
        negative_fitness.push_back(negative_subpopulation_fitness);
    }

    return negative_fitness;
}

void GeneticPatternExtractor::select_new_negative_population(
        std::vector<std::vector<std::vector<Edge *>>> &population,
        std::vector<std::vector<double>> &fitness,
        std::vector<int> &subpopulation_widths,
        std::vector<std::vector<std::vector<double>>> &negative_partial_fitnesses,
        std::vector<std::vector<std::vector<int>>> &negative_partial_counts)
{
    std::vector<std::vector<std::vector<Edge *>>> new_population;
    std::vector<std::vector<std::vector<double>>> new_partial_fitnesses;
    std::vector<std::vector<std::vector<int>>> new_partial_counts;

    for (int subpopulation_index = 0; subpopulation_index < population.size(); ++subpopulation_index)
    {

        std::vector<std::vector<Edge *>> new_subpopulation;
        std::vector<std::vector<double>> new_subpopulation_partial_fitness;
        std::vector<std::vector<int>> new_subpopulation_partial_counts;

        if (population.at(subpopulation_index).size() > 0)
        {
            double min_fitness = fitness.at(subpopulation_index).at(0);
            for (int i = 0; i < fitness.at(subpopulation_index).size(); ++i) {
                if (fitness.at(subpopulation_index).at(i) < min_fitness)
                {
                    min_fitness = fitness.at(subpopulation_index).at(i);
                }
            }
            std::vector<int> selected_individuals_indices;
            if (min_fitness <= 0.0)
            {
                // if there are some negative fitness values, shift all values so that minimum is zero + small smoothing
                std::vector<double> weights;
                for (auto& fit : fitness.at(subpopulation_index))
                {
                    weights.push_back(fit - min_fitness + 0.001);
                }
                selected_individuals_indices = random_generator->stochastic_universal_sampling(
                        subpopulation_widths.at(subpopulation_index), weights);
            }
            else
            {
                // otherwise just use the fitness values as they are
                selected_individuals_indices = random_generator->stochastic_universal_sampling(
                        subpopulation_widths.at(subpopulation_index), fitness.at(subpopulation_index));
            }

            for (auto& ind : selected_individuals_indices)
            {
                new_subpopulation.push_back(population.at(subpopulation_index).at(ind));
                new_subpopulation_partial_fitness.push_back(negative_partial_fitnesses.at(subpopulation_index).at(ind));
                new_subpopulation_partial_counts.push_back(negative_partial_counts.at(subpopulation_index).at(ind));
            }
        }

        new_population.push_back(new_subpopulation);
        new_partial_fitnesses.push_back(new_subpopulation_partial_fitness);
        new_partial_counts.push_back(new_subpopulation_partial_counts);
    }

    // replace old vectors with the new ones:
    population = new_population;
    negative_partial_fitnesses = new_partial_fitnesses;
    negative_partial_counts = new_partial_counts;
}


std::vector<double> GeneticPatternExtractor::select_best_negative_individual_for_each_positive(
        std::vector<std::vector<std::vector<Edge *>>> &negative_population,
        std::vector<std::vector<double>> &negative_fitness)
{
    std::vector<double> negative_individuals_score;
    for (int subpopulation_index = 0; subpopulation_index < negative_population.size(); ++subpopulation_index)
    {
        double best_score = 0.0;
        for (int individual_index = 0; individual_index < negative_population.at(subpopulation_index).size(); ++individual_index)
        {
            if (negative_fitness.at(subpopulation_index).at(individual_index) > best_score)
            {
                best_score = negative_fitness.at(subpopulation_index).at(individual_index);
            }
        }
        negative_individuals_score.push_back(best_score);
    }
    return negative_individuals_score;
}

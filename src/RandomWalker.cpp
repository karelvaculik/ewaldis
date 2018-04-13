//
// Created by karbal on 29.12.17.
//

#include <iostream>
#include "RandomWalker.h"
#include "commonutilstemplated.h"

RandomWalker::RandomWalker(bool use_vertex_attributes,
                           timestamp_t time_unit_primary,
                           timestamp_t time_unit_secondary,
                           int random_walks,
                           double prob_restart,
                           RandomGenerator * random_generator) : use_vertex_attributes(use_vertex_attributes),
                                                               time_unit_primary(time_unit_primary),
                                                               time_unit_secondary(time_unit_secondary),
                                                               random_walks(random_walks),
                                                               prob_restart(prob_restart),
                                                               random_generator(random_generator)
{
}


bool RandomWalker::prob_succeed(double prob)
{
    double value = rand() / (double) RAND_MAX;
    return value < prob;
}


double RandomWalker::expon_pdf(double x, double lamb, timestamp_t time_unit)
{
    if (x < 0) {
        return 0.0;
    } else {
        return lamb * exp(-lamb * x / time_unit);
    }
}

double RandomWalker::edge_similarity(DynamicGraph * graph, Edge * primary_edge, Edge * secondary_edge,
                                    timestamp_t secondary_event_start_time,
                                    timestamp_t secondary_edge_expected_timestamp)
{
    if (primary_edge->get_direction() != secondary_edge->get_direction() ||
        secondary_edge->get_timestamp() > secondary_event_start_time) {
        return 0.0;
    }

    double attribute_distance = 0.0;

    // compute the mixed-euclidean distance between the attributes of the edges, first just the differences
    for (auto const &attribute_name_type_pair : graph->get_edge_schema()) {
        if (attribute_name_type_pair.second == AttributeType::NOMINAL) {
            // if the attribute is nominal, just compute the 0/1 distance
            int f = (int) primary_edge->get_attributes().at(attribute_name_type_pair.first);
            int s = (int) secondary_edge->get_attributes().at(attribute_name_type_pair.first);
            attribute_distance += double(f != s);
        } else if (attribute_name_type_pair.second == AttributeType::NUMERIC) {
            // if the attribute is numeric, it is the difference squared
            attribute_distance += pow(primary_edge->get_attributes().at(attribute_name_type_pair.first) -
                                      secondary_edge->get_attributes().at(attribute_name_type_pair.first), 2);
        }
    }
    // now also for vertices
    if (this->use_vertex_attributes) {
        for (auto const &attribute_name_type_pair : graph->get_vertex_schema()) {
            if (attribute_name_type_pair.second == AttributeType::NOMINAL) {
                // if the attribute is nominal, just compute the 0/1 distance

                int from_f = (int) graph->get_vertices().at(primary_edge->get_from_vertex_id())->get_attributes().at(
                        attribute_name_type_pair.first);
                int from_s = (int) graph->get_vertices().at(secondary_edge->get_from_vertex_id())->get_attributes().at(
                        attribute_name_type_pair.first);
                attribute_distance += double(from_f != from_s);
                int to_f = (int) graph->get_vertices().at(primary_edge->get_to_vertex_id())->get_attributes().at(
                        attribute_name_type_pair.first);
                int to_s = (int) graph->get_vertices().at(secondary_edge->get_to_vertex_id())->get_attributes().at(
                        attribute_name_type_pair.first);
                attribute_distance += double(to_f != to_s);
            } else if (attribute_name_type_pair.second == AttributeType::NUMERIC) {
                // if the attribute is numeric, it is the difference squared
                attribute_distance += pow(
                        graph->get_vertices().at(primary_edge->get_from_vertex_id())->get_attributes().at(
                                attribute_name_type_pair.first)
                        - graph->get_vertices().at(secondary_edge->get_from_vertex_id())->get_attributes().at(
                                attribute_name_type_pair.first), 2);
                attribute_distance += pow(graph->get_vertices().at(primary_edge->get_to_vertex_id())->get_attributes().at(
                        attribute_name_type_pair.first)
                                          - graph->get_vertices().at(
                        secondary_edge->get_to_vertex_id())->get_attributes().at(attribute_name_type_pair.first), 2);
            }
        }
    }

    attribute_distance = sqrt(attribute_distance);
    // compute maximum attribute distance so that we can normalize the attribute distance
    double max_attribute_distance;
    if (use_vertex_attributes) {
        max_attribute_distance = sqrt(graph->get_edge_schema().size() + 2.0 * graph->get_vertex_schema().size());
    } else {
        max_attribute_distance = sqrt(graph->get_edge_schema().size());
    }
    // compute attribute similarity and time similarity and then the total similarity from those two
    double attribute_similarity = std::max(0.1, 1.0 - attribute_distance / max_attribute_distance);

    double time_similarity = expon_pdf(abs(secondary_edge_expected_timestamp - secondary_edge->get_timestamp()),
                                      1.0, this->time_unit_secondary);


    double total_similarity = 2.0 * attribute_similarity * time_similarity / (attribute_similarity + time_similarity);

    return total_similarity;
}

std::pair<void *, double>
RandomWalker::select_primary_edge(DynamicGraph * graph, int current_node_id, timestamp_t start_time,
                                  std::set<int> already_visited_edges)
{

    // first, take only edges that haven't been visited yet:
    std::vector<Edge *> allowedEdges;
    std::vector<double> probs;

    for (int j = 0; j < graph->get_adjacency_list().at(current_node_id).size(); ++j)
    {
        Edge * edge = graph->get_adjacency_list().at(current_node_id).at(j);
        if (already_visited_edges.find(edge->get_original_edge_id()) == already_visited_edges.end()) {
            // if the edge is not in already_visited_edges, use it
            double prob = expon_pdf(start_time - edge->get_timestamp(), 1.0, this->time_unit_primary);
            if (prob > 0.0) {
                allowedEdges.push_back(edge);
                probs.push_back(prob);
            }

        }
    }

    if (allowedEdges.size() == 0) {
        return std::make_pair(nullptr, 0.0);
    }
    int selected_index = random_generator->generate_random_int(0, allowedEdges.size());

    return std::make_pair(allowedEdges[selected_index], probs[selected_index]);
}

std::pair<void *, double>
RandomWalker::select_secondary_edge(DynamicGraph * graph, int current_node_id, timestamp_t secondary_event_start_time,
                                    timestamp_t secondary_edge_expected_timestamp, Edge * primary_edge,
                                    std::set<int> already_visited_edges)
{
    std::vector<Edge *> allowedEdges;
    std::vector<double> probs;

    for (int j = 0; j < graph->get_adjacency_list().at(current_node_id).size(); ++j)
    {
        Edge * edge = graph->get_adjacency_list().at(current_node_id).at(j);
        if (already_visited_edges.find(edge->get_original_edge_id()) == already_visited_edges.end()) {
            // if the edge is not in already_visited_edges, use it
            double prob = edge_similarity(graph, primary_edge, edge, secondary_event_start_time,
                                         secondary_edge_expected_timestamp);
            allowedEdges.push_back(edge);
            probs.push_back(prob);

        }
    }

    double prob_sum = 0.0;
    for (auto &prob : probs) {
        prob_sum += prob;
    }

    if (fabs(prob_sum - 0.0) < 10e-8) {
        return std::make_pair(nullptr, 0.0);
    }

    int selected_index = random_generator->generate_random_int_from_distribution(probs);

    return std::make_pair(allowedEdges[selected_index], probs[selected_index]);
}

void RandomWalker::perform_one_walk(std::vector<DynamicGraph> & graph_instances, std::vector<int> & current_nodes,
                                    std::vector<timestamp_t > & start_nodes_times, int n_positive_instances,
                                    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_positive,
                                    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_negative,
                                    std::map<std::pair<int, int>, double> & edge_priors,
                                    std::vector<std::set<int>> & already_visited_edges)
{
    int n_instances = current_nodes.size();

    // copy the elements from already_visited_edges
    std::vector<std::set<int>> used_edge_sets = already_visited_edges;

    // generate vector of n_instances numbers
    std::vector<int> instances_alive = generate_integer_sequence(n_instances);

    std::vector<int> positive_instances_alive = generate_integer_sequence(n_positive_instances);

    while (!prob_succeed(this->prob_restart)) {
        if (positive_instances_alive.size() <= 1) {
            break;
        }

        int rand_int = random_generator->generate_random_int(0, positive_instances_alive.size());
        int primary_occurrence_index = positive_instances_alive[rand_int];


        std::pair<void *, double> primary_edge_prob = select_primary_edge(&graph_instances[primary_occurrence_index],
                                                                          current_nodes[primary_occurrence_index],
                                                                          start_nodes_times[primary_occurrence_index],
                                                                          used_edge_sets[primary_occurrence_index]);

        if (primary_edge_prob.first == nullptr) {
            break;
        }

        Edge *selected_primary_edge = (Edge *) primary_edge_prob.first;
        double edge_prob = primary_edge_prob.second;

        edge_priors[std::make_pair(primary_occurrence_index, selected_primary_edge->get_edge_id())] = edge_prob;

        timestamp_t primary_edge_timestamp_difference =
                start_nodes_times[primary_occurrence_index] - selected_primary_edge->get_timestamp();

        used_edge_sets[primary_occurrence_index].insert(selected_primary_edge->get_original_edge_id());
        current_nodes[primary_occurrence_index] = selected_primary_edge->get_to_vertex_id();

        std::vector<int> new_instances_alive;
        std::vector<int> new_positive_instances_alive;
        new_instances_alive.push_back(primary_occurrence_index);
        new_positive_instances_alive.push_back(primary_occurrence_index);

        for (auto &instance_index : instances_alive) {

            if (instance_index != primary_occurrence_index) {
                timestamp_t secondary_edge_expected_timestamp =
                        start_nodes_times[instance_index] - primary_edge_timestamp_difference;

                std::pair<void *, double> secondary_edge_prob = select_secondary_edge(&graph_instances[instance_index],
                                                                                      current_nodes[instance_index],
                                                                                      start_nodes_times[instance_index],
                                                                                      secondary_edge_expected_timestamp,
                                                                                      selected_primary_edge,
                                                                                      used_edge_sets[instance_index]);

                if (secondary_edge_prob.first != nullptr) {


                    Edge *selected_secondary_edge = (Edge *) secondary_edge_prob.first;
                    double similarity = secondary_edge_prob.second;

                    std::pair<int, int> primary_pair = std::make_pair(primary_occurrence_index,
                                                                      selected_primary_edge->get_edge_id());
                    std::pair<int, int> secondary_pair = std::make_pair(instance_index,
                                                                        selected_secondary_edge->get_edge_id());
                    //  we just save the similarity (don't accumulate it)
                    if (instance_index < n_positive_instances)
                    {
                        // save similarity for positive instance (NEW: use smaller index as primary key)
                        if (primary_occurrence_index < instance_index) {
                            if (edge_pairs_dictionary_positive.find(primary_pair) ==
                                edge_pairs_dictionary_positive.end()) {
                                edge_pairs_dictionary_positive[primary_pair] = std::map<std::pair<int, int>, double>();
                            }
                            edge_pairs_dictionary_positive.at(primary_pair)[secondary_pair] = similarity;
                        } else {
                            if (edge_pairs_dictionary_positive.find(secondary_pair) ==
                                edge_pairs_dictionary_positive.end()) {
                                edge_pairs_dictionary_positive[secondary_pair] = std::map<std::pair<int, int>, double>();
                            }
                            edge_pairs_dictionary_positive.at(secondary_pair)[primary_pair] = similarity;

                        }
                    }
                    else
                    {
                        if (edge_pairs_dictionary_negative.find(primary_pair) == edge_pairs_dictionary_negative.end()) {
                            edge_pairs_dictionary_negative[primary_pair] = std::map<std::pair<int, int>, double>();
                        }
                        edge_pairs_dictionary_negative.at(primary_pair)[secondary_pair] = similarity;
                    }

                    current_nodes[instance_index] = selected_secondary_edge->get_to_vertex_id();
                    new_instances_alive.push_back(instance_index);
                    if (instance_index < n_positive_instances) {
                        new_positive_instances_alive.push_back(instance_index);
                    }
                    used_edge_sets[instance_index].insert(selected_secondary_edge->get_original_edge_id());
                }
            }
        }
        instances_alive = new_instances_alive;
        positive_instances_alive = new_positive_instances_alive;

    }
}


Suitabilities RandomWalker::prepare_scores(
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_positive,
        std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_negative,
        std::map<std::pair<int, int>, double> & edge_priors)
{
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> new_edge_pairs_dictionary_positive;
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> new_edge_pairs_dictionary_negative;

    // process positive
    for (auto const &pos_pair : edge_pairs_dictionary_positive) {

        std::pair<int, int> first = pos_pair.first;


        if (edge_priors.find(first) == edge_priors.end()) {
            // TODO: in one_walk store also the opposite edge in priors and remove this part of code
            continue;
        }

        double A = edge_priors.at(first);
        for (auto const &pos_pair_inner : pos_pair.second) {
            std::pair<int, int> second = pos_pair_inner.first;
            double B;
            if (edge_priors.find(second) != edge_priors.end()) {
                B = edge_priors.at(second);
            } else {
                B = A;
            }
            double C = pos_pair_inner.second;

            double score = (A + B) * C / 2.0;
            if (score < 0.1) {
                // if the score is too low, just don't use it
                continue;
            }

            if (new_edge_pairs_dictionary_positive.find(first) == new_edge_pairs_dictionary_positive.end()) {
                new_edge_pairs_dictionary_positive[first] = std::map<std::pair<int, int>, double>();
            }
            new_edge_pairs_dictionary_positive.at(first)[second] = score;

            if (new_edge_pairs_dictionary_positive.find(second) == new_edge_pairs_dictionary_positive.end()) {
                new_edge_pairs_dictionary_positive[second] = std::map<std::pair<int, int>, double>();
            }
            new_edge_pairs_dictionary_positive.at(second)[first] = score;

        }
    }
    // process negative
    for (auto const &neg_pair : edge_pairs_dictionary_negative) {

        std::pair<int, int> first = neg_pair.first;
        double A = edge_priors.at(first);
        for (auto const &neg_pair_inner : neg_pair.second) {
            std::pair<int, int> second = neg_pair_inner.first;
            double C = neg_pair_inner.second;
            double score = A * C;
            if (score < 0.1) {
                // if the score is too low, just don't use it
                continue;
            }

            if (new_edge_pairs_dictionary_negative.find(first) == new_edge_pairs_dictionary_negative.end()) {
                new_edge_pairs_dictionary_negative[first] = std::map<std::pair<int, int>, double>();
            }
            new_edge_pairs_dictionary_negative.at(first)[second] = score;
        }
    }

    return Suitabilities(new_edge_pairs_dictionary_positive, new_edge_pairs_dictionary_negative);
}


Suitabilities RandomWalker::compute_suitabilities(std::vector<DynamicGraph> & graph_instances,
                                                  std::vector<std::vector<int>> &positive_event_vertices,
                                                  std::vector<timestamp_t> &positive_event_times,
                                                  std::vector<std::vector<int>> &positive_event_edges,
                                                  std::vector<std::vector<int>> &negative_event_vertices,
                                                  std::vector<timestamp_t> &negative_event_times,
                                                  std::vector<std::vector<int>> &negative_event_edges)
{
    // positive_event_vertices for each variant, there is a list of vertex ids (one for each instance)
    // same for positive_event_edges

    // how many positive instances there are:
    int n_positive_instances = positive_event_vertices.at(0).size();

    // for each instance, we have a set of edges that are currently occupied
    // initially, there are edges of the initial patterns
    std::vector<std::set<int>> already_visited_edges;

    // init the vector of sets:
    for (int l = 0; l < graph_instances.size(); ++l)
    {
        already_visited_edges.push_back(std::set<int>());
    }

    // process the positive event edges
    for (int i = 0; i < positive_event_edges.size(); ++i) {
        // for each variant
        for (int j = 0; j < positive_event_edges.at(i).size(); ++j) {
            // for each instance
            already_visited_edges.at(j).insert(positive_event_edges.at(i).at(j));
        }
    }
    // process the negative event edges
    for (int i = 0; i < negative_event_edges.size(); ++i) {
        // for each variant
        for (int j = 0; j < negative_event_edges.at(i).size(); ++j) {
            // for each instance
            already_visited_edges.at(j + n_positive_instances).insert(negative_event_edges.at(i).at(j));
        }
    }

    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_positive;
    std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> edge_pairs_dictionary_negative;
    std::map<std::pair<int, int>, double> edge_priors;

    int which_nodes;

    std::vector<timestamp_t > start_times;
    start_times.insert(start_times.end(), positive_event_times.begin(), positive_event_times.end());
    start_times.insert(start_times.end(), negative_event_times.begin(), negative_event_times.end());

    for (int k = 0; k < this->random_walks; ++k) {
        which_nodes = random_generator->generate_random_int(0, positive_event_vertices.size());

        std::vector<int> start_nodes;
        start_nodes.insert(start_nodes.end(), positive_event_vertices.at(which_nodes).begin(),
                           positive_event_vertices.at(which_nodes).end());
        start_nodes.insert(start_nodes.end(), negative_event_vertices.at(which_nodes).begin(),
                           negative_event_vertices.at(which_nodes).end());

        perform_one_walk(graph_instances, start_nodes, start_times, n_positive_instances, edge_pairs_dictionary_positive,
                         edge_pairs_dictionary_negative, edge_priors, already_visited_edges);
    }

    // finally, compute the suitabilities from the dictionaries
    return prepare_scores(edge_pairs_dictionary_positive, edge_pairs_dictionary_negative, edge_priors);
}
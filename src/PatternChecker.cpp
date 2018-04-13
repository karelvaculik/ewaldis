//
// Created by karel on 14.1.18.
//

#include <iostream>
#include "PatternChecker.h"
#include "RandomWalker.h"

#include "commonutilstemplated.h"



PatternChecker::PatternChecker(bool use_vertex_attributes,
               timestamp_t time_unit_secondary,
               int random_walks,
               RandomGenerator * random_generator) :
    use_vertex_attributes(use_vertex_attributes),
    time_unit_secondary(time_unit_secondary),
    random_walks(random_walks),
    random_generator(random_generator)
{
}

double PatternChecker::check_pattern_in_instance(DynamicGraph * graph, Pattern * pattern,
                                 std::vector<int> & graph_starting_vertices,
                                 timestamp_t graph_starting_timestamp,
                                 std::vector<double> & edge_weights)
{

    if (pattern->get_pattern_vertex_pairs().size() == 0) {
        return 1.0;
    }
    int number_of_instances_in_pattern = (int) pattern->get_attributes().at(0).size();

    // prepare adjacency list from pattern edges
    std::map<int, std::vector<int>> adjacency_list;
    // we also remember which edge is linked to each such adjacency information - it is necessary for score computation
    std::map<int, std::vector<int>> adjacency_links_to_edges;
    // whether this is the original direction or the added (opposite) one:
    std::map<int, std::vector<bool>> adjacency_list_is_opposite;
    for (int i = 0; i < pattern->get_pattern_vertex_pairs().size(); ++i)
    {
        int src = pattern->get_pattern_vertex_pairs().at(i).first;
        int dst = pattern->get_pattern_vertex_pairs().at(i).second;
        if (adjacency_list.find(src) == adjacency_list.end())
        {
            adjacency_list[src] = std::vector<int>();
            adjacency_links_to_edges[src] = std::vector<int>();
            adjacency_list_is_opposite[src] = std::vector<bool>();
        }
        adjacency_list.at(src).push_back(dst);
        adjacency_links_to_edges.at(src).push_back(i);
        adjacency_list_is_opposite.at(src).push_back(false);

        if (adjacency_list.find(dst) == adjacency_list.end())
        {
            adjacency_list[dst] = std::vector<int>();
            adjacency_links_to_edges[dst] = std::vector<int>();
            adjacency_list_is_opposite[dst] = std::vector<bool>();
        }
        adjacency_list.at(dst).push_back(src);
        adjacency_links_to_edges.at(dst).push_back(i);
        adjacency_list_is_opposite.at(dst).push_back(true);
    }

    // not all starting vertices occur in the pattern, so use 0 selection probability for those missing
    std::vector<bool> starting_vertices_usable(graph_starting_vertices.size());
    std::fill(starting_vertices_usable.begin(), starting_vertices_usable.end(), false);
    for (auto& pattern_vertex_pair : pattern->get_pattern_vertex_pairs())
    {
        if (pattern_vertex_pair.first < graph_starting_vertices.size())
        {
            starting_vertices_usable[pattern_vertex_pair.first] = true;
        }
        if (pattern_vertex_pair.second < graph_starting_vertices.size())
        {
            starting_vertices_usable[pattern_vertex_pair.second] = true;
        }
    }

    // use this simple random walker in order to able to compute the similarities
    RandomWalker random_walker = RandomWalker(use_vertex_attributes, 0.0, time_unit_secondary, 0, 0.0, random_generator);



    std::vector<double> all_total_scores;

    for (int pattern_instance_index = 0; pattern_instance_index < number_of_instances_in_pattern; ++pattern_instance_index)
    {
        std::vector<double> walks_scores;
        for (int i = 0; i < random_walks; ++i) {
            // try one mapping of the pattern

            // what is the score of this walk and how many steps we did
            double walk_score = 0.0;
            int walk_length = 0;

            // we don't allow to use one pattern edge many times (both in pattern and graph) in one random walk
            std::set<int> already_visited_pattern_edges;
            // here we keep the original ids used in the graph in this single random walk
            std::set<int> already_visited_graph_edges;

            // here we keep the vertices from which we can go one edge
            std::vector<int> available_vertices;

            // here we keep the vertices that have been already added to available_vertices
            // so that we do not add them there again (once they are removed)
            std::set<int> already_tried_vertices;

            // keep the mapping from pattern vertices to graph vertices and use it to check the vertex consistency
            std::map<int, int> vertex_mapping;
            for (int j = 0; j < graph_starting_vertices.size(); ++j)
            {
                if (starting_vertices_usable.at(j))
                {
                    available_vertices.push_back(j);
                    already_tried_vertices.insert(j);
//                    vertex_mapping.emplace(j, graph_starting_vertices.at(j));
                    vertex_mapping[j] = graph_starting_vertices.at(j);
                }
            }

            // now perform the occupation of the graph by the pattern
            while (available_vertices.size() > 0)
            {
                // while we have some available vertices (i.e. available edges)
                // pick one edge randomly
                int current_pattern_vertex_index = random_generator->generate_random_int(0, available_vertices.size());
                int current_pattern_vertex = available_vertices.at(current_pattern_vertex_index);

                int current_graph_vertex = vertex_mapping.at(current_pattern_vertex);

                std::vector<int> which_adjacent_edges_could_be_used;
                for (int j = 0; j < adjacency_links_to_edges.at(current_pattern_vertex).size(); ++j)
                {
                    int considered_element = adjacency_links_to_edges.at(current_pattern_vertex).at(j);
                    if (already_visited_pattern_edges.find(considered_element) == already_visited_pattern_edges.end())
                    {
                        which_adjacent_edges_could_be_used.push_back(j);
                    }
                }

                if (which_adjacent_edges_could_be_used.size() == 0)
                {
                    // if there are no available pattern edges going from this vertex,
                    // remove it from our set and try another round
                    available_vertices.erase(available_vertices.begin() + current_pattern_vertex_index);
                    continue;
                }

                // if there are some usable edges, select one at random
                int index_in_adjacency_list_j = random_generator->generate_random_int(0, which_adjacent_edges_could_be_used.size());
                int index_in_adjacency_list = which_adjacent_edges_could_be_used.at(index_in_adjacency_list_j);

                int pattern_edge_index = adjacency_links_to_edges.at(current_pattern_vertex).at(index_in_adjacency_list);

                // increase the counts
                if (pattern->get_attributes().at(pattern_edge_index).at(pattern_instance_index).size() > 0)
                {
                    // don't count walk length if there nothing here (size == 0)
                    ++walk_length;
                }
                already_visited_pattern_edges.insert(pattern_edge_index);

                int pattern_to_vertex = adjacency_list.at(current_pattern_vertex).at(index_in_adjacency_list);

                if (already_tried_vertices.find(pattern_to_vertex) == already_tried_vertices.end())
                {
                    // if the vertex is not in the available_vertices, add it there
                    already_tried_vertices.insert(pattern_to_vertex);
                }

                std::vector<Edge *> graph_allowed_edges;
                bool pattern_to_vertex_in_vertex_mapping = vertex_mapping.find(pattern_to_vertex) != vertex_mapping.end();
                for (auto& e : graph->get_adjacency_list().at(current_graph_vertex))
                {
                    if (already_visited_graph_edges.find(e->get_original_edge_id()) == already_visited_graph_edges.end())
                    {
                        if (pattern_to_vertex_in_vertex_mapping)
                        {
                            // if patter_to_vertex is already in vertex_mapping, check also consistency:
                            // remove edges that have inconsistent vertices (we check dst vertex)
                            if (e->get_to_vertex_id() == vertex_mapping.at(pattern_to_vertex))
                            {
                                graph_allowed_edges.push_back(e);
                            }
                        }
                        else
                        {
                            graph_allowed_edges.push_back(e);
                        }
                    }
                }
                if (graph_allowed_edges.size() == 0)
                {
                    // update available vertices and continue by selecting another edge
                    available_vertices.erase(available_vertices.begin() + current_pattern_vertex_index);
                    continue;
                }
                else
                {
                    // select the edge in the graph and move forward both in the graph and pattern
                    // only attributes and direction are necessary for edge similarity function

                    std::vector<double> probs;
                    if (pattern->get_attributes().at(pattern_edge_index).at(pattern_instance_index).size() > 0)
                    {
                        // don't count walk length if there None here
                        Edge primary_edge = Edge(0, 0, 0, 0.0, pattern->get_attributes().at(pattern_edge_index).at(pattern_instance_index),
                                    pattern->get_directions().at(pattern_edge_index), constants::INVALID_ID);
                        // now check whether we should take the opposite
                        if (adjacency_list_is_opposite.at(current_pattern_vertex).at(index_in_adjacency_list))
                        {
                            // take the opposite
                            Edge opposite_edge = primary_edge.create_opposite_edge(graph->is_undirected(), 0);
                            primary_edge = opposite_edge;
                        }
                        for(auto& secondary_edge : graph_allowed_edges)
                        {
                            probs.push_back(random_walker.edge_similarity(graph, &primary_edge, secondary_edge,
                                                          graph_starting_timestamp,
                                                          graph_starting_timestamp - pattern->get_timestamps().at(pattern_edge_index).at(pattern_instance_index)));
                        }
                    }
                    else
                    {
                        // if the pattern edge is None, you can use any edge from the graph
                        std::fill(probs.begin(), probs.end(), 1.0);
                    }

                    double probs_sum = 0.0;
                    for (auto& pr : probs)
                    {
                        probs_sum += pr;
                    }

                    if (probs_sum == 0.0)
                    {
                        // we cannot select anything, so continue by selecting another edge
                        // also update the available vertices
                        available_vertices.erase(available_vertices.begin() + current_pattern_vertex_index);
                        continue;
                    }
                    else
                    {
                        // we are able to select an edge in the graph, so select one and compute the similarity
                        int selected_graph_edge_index = random_generator->generate_random_int_from_distribution(probs);
                        Edge * selected_graph_edge = graph_allowed_edges.at(selected_graph_edge_index);


                        already_visited_graph_edges.insert(selected_graph_edge->get_original_edge_id());

                        if (pattern->get_attributes().at(pattern_edge_index).at(pattern_instance_index).size() > 0)
                        {
                            // update score only if the edge wasn't None
                            walk_score += probs.at(selected_graph_edge_index);
                        }
                    }
                }
            }

            if (walk_length > 0)
            {
                walks_scores.push_back(walk_score / walk_length);
            }
            else
            {
                walks_scores.push_back(0.0);
            }
        }
        // now take the maximum walk score and save it as a score of this pattern instance:
        double max_score = walks_scores.at(0);
        for (auto& sc : walks_scores)
        {
            if (sc > max_score)
            {
                max_score = sc;
            }
        }
        all_total_scores.push_back(max_score);
    }

    // return mean of all_total_scores
    double sum_of_all = 0.0;
    for (auto& sc : all_total_scores)
    {
        sum_of_all += sc;
    }
    return sum_of_all / all_total_scores.size();
}
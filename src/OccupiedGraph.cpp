//
// Created by karbal on 2.1.18.
//

#include <iostream>

#include "OccupiedGraph.h"

#include "commonutilstemplated.h"


using namespace std;


OccupiedGraph::~OccupiedGraph() {}


OccupiedGraph::OccupiedGraph(std::vector<DynamicGraph> & graph_instances,
                             const std::vector<std::vector<int>> & positive_event_vertices,
                             const std::vector<std::vector<int>> & positive_event_edges,
                             const std::map<std::pair<int, int>, std::map<std::pair<int, int>, double>> & edge_pairs_dictionary_positive) :
        original_graph_instances(graph_instances),
        edge_pairs_dictionary_positive(edge_pairs_dictionary_positive)
{
    // here we keep all (INSTANCE, ORIGINAL EDGE ID) pairs of already selected edges, so we can check that
    // we do not select an edge twice for an instance
    // for each edge combination:
    for (int i = 0; i < positive_event_edges.size(); ++i)
    {
        // for each instance in the edge combination
        for (int j = 0; j < positive_event_edges.at(i).size(); ++j)
        {
            // j = instance_id
            all_occupied_original_edges_ids_with_instances.insert(std::make_pair(j, positive_event_edges.at(i).at(j)));
        }
    }

    // how many positive instances there:
    n_positive = positive_event_vertices.at(0).size();
    vertex_sets_current_id = 0;
    // add the initial pattern nodes to vertex sets mappings:
    // for each vertex combination:
    for (int k = 0; k < positive_event_vertices.size(); ++k)
    {
        // for each instance in the vertex combination
        for (int i = 0; i < positive_event_vertices.at(k).size(); ++i)
        {
            // i = instance_id
            if (vertex_mapping.find( i ) == vertex_mapping.end())
            {
                // initialize the maps:
                vertex_mapping[i] = std::map<int, int>();
                vertex_mapping_reversed[i] = std::map<int, int>();
            }
            vertex_mapping.at(i)[positive_event_vertices.at(k).at(i)] = vertex_sets_current_id;
            vertex_mapping_reversed.at(i)[vertex_sets_current_id] = positive_event_vertices.at(k).at(i);
        }
        ++vertex_sets_current_id;
    }
    prepare_initial_edge_combination(combinations, to_vertex_id_keys);
}


void OccupiedGraph::prepare_initial_edge_combination(std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & combinations,
                                                     std::map<std::pair<int, int>, std::set<std::pair<int, int>>> & to_vertex_id_keys)
{
    // get all vertex set ids, so that we can process the nodes in the same order across all instances
    std::vector<int> vertex_set_ids;
    for (auto const& element : vertex_mapping_reversed.at(0)) {
        vertex_set_ids.push_back(element.first);
    }
    for (int instance_id = 0; instance_id < n_positive; ++instance_id)
    {
        // for each instance, prepare edges by using the starting vertices
        for (auto const& from_vertex_set_id : vertex_set_ids)
        {
            int from_vertex_id = vertex_mapping_reversed.at(instance_id).at(from_vertex_set_id);
            // process all edges going from this vertex:
            std::vector<Edge *> adjacent_edges = original_graph_instances.at(instance_id).get_adjacency_list().at(from_vertex_id);
            for (Edge * adjacent_edge : adjacent_edges)
            {
                int original_edge_id = adjacent_edge->get_original_edge_id();
                int edge_id = adjacent_edge->get_edge_id();
                std::pair<int, int> original_id_pair = std::make_pair(instance_id, original_edge_id);
                std::pair<int, int> id_pair = std::make_pair(instance_id, edge_id);
                bool is_occupied = (all_occupied_original_edges_ids_with_instances.find(original_id_pair) != all_occupied_original_edges_ids_with_instances.end());
                bool is_candidate = (original_edges_ids_with_instances_in_combinations.find(original_id_pair) != original_edges_ids_with_instances_in_combinations.end());
                bool is_in_dictionary = (edge_pairs_dictionary_positive.find(id_pair) != edge_pairs_dictionary_positive.end());
                if (is_occupied || is_candidate || !is_in_dictionary)
                {
                    // if the edge is already occupied or in combinations or it is not in the dictionary, skip it
                    continue;
                }
                else
                {
                    // otherwise add it to the combinations
                    int to_vertex_id = adjacent_edge->get_to_vertex_id();
                    int to_vertex_set_id = get_vertex_set(instance_id, to_vertex_id);
                    std::pair<int, int> to_vertex_pair = std::make_pair(instance_id, to_vertex_id);
                    if (to_vertex_id_keys.find(to_vertex_pair) == to_vertex_id_keys.end())
                    {
                        to_vertex_id_keys[to_vertex_pair] = std::set<std::pair<int, int>>();
                    }
                    std::pair<int, int> from_to_id_pair = std::make_pair(from_vertex_set_id, to_vertex_set_id);
                    to_vertex_id_keys.at(to_vertex_pair).insert(from_to_id_pair);

                    add_to_combinations(combinations, from_to_id_pair, instance_id, adjacent_edge, n_positive);

                    original_edges_ids_with_instances_in_combinations.insert(original_id_pair);

                }
            }

        }
    }
}


int OccupiedGraph::get_vertex_set(int instance_id, int vertex_id)
{

    // Returns the vertex_set id of a given vertex in a given instance.
    // If the vertex is not mapped yet, it returns INVALID_ID.
    if (vertex_mapping.find(instance_id) == vertex_mapping.end())
    {
        throw std::runtime_error("Instance ID is not found in the vertex_mapping");
    }
    else
    {
        if (vertex_mapping.at(instance_id).find(vertex_id) == vertex_mapping.at(instance_id).end())
        {
            return constants::INVALID_ID;
        }
        else
        {
            return vertex_mapping.at(instance_id).at(vertex_id);
        }
    }
}

void OccupiedGraph::add_to_combinations(std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & combinations,
                                        std::pair<int, int> from_to_id_pair,
                                        int instance_id, Edge * edge, int n_positive)
{
    // from_to_id_pair = pair of IDS (from-vertex, to-vertex)
    if (combinations.find(from_to_id_pair) == combinations.end())
    {
        // create an empty vector of edges for each instance_id and use this vector of empty vectors as a new combination
        std::vector<std::vector<Edge *>> empty_combination;
        for (int i = 0; i < n_positive; ++i) {
            empty_combination.push_back(std::vector<Edge *>());
        }
        combinations[from_to_id_pair] = empty_combination;
    }
    // now add the edge to vector under the right pair of from/to-vertex ids and instance id:
    combinations.at(from_to_id_pair).at(instance_id).push_back(edge);
}

const std::map<std::pair<int, int>, std::vector<std::vector<Edge *>>> & OccupiedGraph::get_edge_combinations()
{
    return combinations;
}


void OccupiedGraph::occupy_edges(const std::vector<Edge*> & edges_to_be_occupied)
{
    // Occupies a set of edges and updates the edge combinations.
    // edges_to_be_occupied: map from instance_id to Edge (because there may not be edges for all instances


    // occupy vertices and edges
    // use the first edge that is not None to get the vertex sets
    int from_vertex_set_id = constants::INVALID_ID;
    int to_vertex_set_id = constants::INVALID_ID;

    for (int instance_id = 0; instance_id < edges_to_be_occupied.size(); ++instance_id) {
        if (edges_to_be_occupied.at(instance_id) != nullptr)
        {
            from_vertex_set_id = get_vertex_set(instance_id, edges_to_be_occupied.at(instance_id)->get_from_vertex_id());
            to_vertex_set_id = get_vertex_set(instance_id, edges_to_be_occupied.at(instance_id)->get_to_vertex_id());
            break;
        }
    }

    for (int instance_id = 0; instance_id < edges_to_be_occupied.size(); ++instance_id)
    {
        if (edges_to_be_occupied.at(instance_id) != nullptr)
        {
            // occupy only the positive edges and only those remove from the possible combinations
            std::pair<int, int> instance_original_id_pair = std::make_pair(instance_id, edges_to_be_occupied.at(instance_id)->get_original_edge_id());
            all_occupied_original_edges_ids_with_instances.insert(instance_original_id_pair);
            original_edges_ids_with_instances_in_combinations.erase(instance_original_id_pair);
            if (to_vertex_set_id == constants::INVALID_ID)
            {
                // if to-vertex == INVALID_ID, it means it is not yet in the vertex mapping, so add it there
                vertex_mapping.at(instance_id)[edges_to_be_occupied.at(instance_id)->get_to_vertex_id()] = vertex_sets_current_id;
                vertex_mapping_reversed.at(instance_id)[vertex_sets_current_id] = edges_to_be_occupied.at(instance_id)->get_to_vertex_id();
            }
        }
    }

    // update the combinations
    // remove the occupied edges
    // here we compute the new combinations for the used vertex-sets pair
    // if the new sub-combination contains some empty set for any instance,
    // we will remove this vertex-sets completely
    std::vector<std::vector<Edge *>> new_sub_combinations;
    std::pair<int, int> from_to_id_pair = std::make_pair(from_vertex_set_id, to_vertex_set_id);

    for (int instance_id = 0; instance_id < combinations.at(from_to_id_pair).size(); ++instance_id) {
        std::vector<Edge *> new_edges_sets;
        int occupied_original_edge_id;
        if (edges_to_be_occupied.at(instance_id) == nullptr)
        {
            occupied_original_edge_id = constants::INVALID_ID;
        }
        else
        {
            occupied_original_edge_id = edges_to_be_occupied.at(instance_id)->get_original_edge_id();
        }
        for (Edge * x : combinations.at(from_to_id_pair).at(instance_id))
        {
            if (x->get_original_edge_id() != occupied_original_edge_id)
            {
                new_edges_sets.push_back(x);
            }
        }
        new_sub_combinations.push_back(new_edges_sets);
    }
    combinations[from_to_id_pair] = new_sub_combinations;

    // add new edges to combinations (in case that a new vertex was added)
    if (to_vertex_set_id == constants::INVALID_ID)
    {
        for (int instance_id = 0; instance_id < edges_to_be_occupied.size(); ++instance_id)
        {
            if (edges_to_be_occupied.at(instance_id) == nullptr)
            {
                continue;
            }
            // to-vertex of the newly occupied edge is a "from-vertex" of edges that are added to combinations
            int to_vertex_id_of_added_edge = edges_to_be_occupied.at(instance_id)->get_to_vertex_id();
            // these adjacent edges will be added as new to combinations
            std::vector<Edge *> adjacent_edges = original_graph_instances.at(instance_id).get_adjacency_list().at(to_vertex_id_of_added_edge);

            // we must also update the key pairs so that they are not (x, None) but (x, y)
            std::set<std::pair<int, int>> new_key_pairs;

            for (Edge *& adjacent_edge : adjacent_edges)
            {
                int original_edge_id = adjacent_edge->get_original_edge_id();
                int edge_id = adjacent_edge->get_edge_id();
                std::pair<int, int> original_id_pair = std::make_pair(instance_id, original_edge_id);
                std::pair<int, int> id_pair = std::make_pair(instance_id, edge_id);

                if (all_occupied_original_edges_ids_with_instances.find(original_id_pair) != all_occupied_original_edges_ids_with_instances.end())
                {
                    // if the edge is in occupied , just skip it
                    continue;
                }
                if (original_edges_ids_with_instances_in_combinations.find(original_id_pair) != original_edges_ids_with_instances_in_combinations.end())
                {
                    // if the edge is already in combinations,
                    // we must change its to-vertex-set from INVALID_ID to an integer
                    // from this perspective, it is from vertex, but from the past perspective,
                    // it was added as to-vertex to INVALID_ID
                    for (auto& key_pair : to_vertex_id_keys.at(std::make_pair(instance_id, adjacent_edge->get_from_vertex_id())))
                    {
                        new_key_pairs.insert(std::make_pair(key_pair.first, vertex_sets_current_id));
                        std::vector<Edge *> new_set_of_edges;
                        for (Edge * set_edge : combinations.at(key_pair).at(instance_id))
                        {
                            if (set_edge->get_to_vertex_id() == adjacent_edge->get_from_vertex_id())
                            {
                                // this edge is to be moved
                                // from_vertex_set stays the same (i.e. key_pair.first),
                                // to_vertex_set is vertex_sets_current_id used earlier
                                std::pair<int, int> from_to_id_pair = std::make_pair(key_pair.first, vertex_sets_current_id);
                                add_to_combinations(combinations, from_to_id_pair, instance_id, set_edge, n_positive);
                            }
                            else
                            {
                                new_set_of_edges.push_back(set_edge);
                            }
                        }
                        // update the set for this pair and instance
                        combinations.at(key_pair).at(instance_id) = new_set_of_edges;
                    }
                }
                else if (edge_pairs_dictionary_positive.find(id_pair) != edge_pairs_dictionary_positive.end())
                {
                    // if the edge is in the dictionary and
                    // the edge is not in the combinations yet, so add it with None as to-vertex-set
                    // there is an exception however and the to-vertex may not be None, but already existing
                    // (in case that the other direction was not added to combinations due to not being in the dictionary)
                    int adjacent_edge_to_vertex_set_id = get_vertex_set(instance_id, adjacent_edge->get_to_vertex_id());

                    std::pair<int, int> instance_to_vertex_pair = std::make_pair(instance_id, adjacent_edge->get_to_vertex_id());
                    if (to_vertex_id_keys.find(instance_to_vertex_pair) == to_vertex_id_keys.end())
                    {
                        to_vertex_id_keys[instance_to_vertex_pair] = std::set<std::pair<int, int>>();
                    }
                    std::pair<int, int> from_to_id_pair = std::make_pair(vertex_sets_current_id, adjacent_edge_to_vertex_set_id);
                    to_vertex_id_keys.at(instance_to_vertex_pair).insert(from_to_id_pair);

                    add_to_combinations(combinations, from_to_id_pair, instance_id, adjacent_edge, n_positive);
                    original_edges_ids_with_instances_in_combinations.insert(original_id_pair);
                }
            }
            to_vertex_id_keys[std::make_pair(instance_id, to_vertex_id_of_added_edge)] = new_key_pairs;
        }
        ++vertex_sets_current_id;
    }

}

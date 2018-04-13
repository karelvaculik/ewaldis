//
// Created by karbal on 29.12.17.
//

#include "Pattern.h"

Pattern::Pattern(std::vector<std::pair<int, int>> pattern_vertex_pairs,
                 std::vector<bool> directions,
                 std::vector<std::vector<timestamp_t>> timestamps,
                 std::vector<std::vector<std::map<std::string, double>>> attributes,
                 std::vector<double> scores) :
        pattern_vertex_pairs(pattern_vertex_pairs),
        directions(directions),
        timestamps(timestamps),
        attributes(attributes),
        scores(scores)
{
}


const std::vector<std::pair<int, int>> & Pattern::get_pattern_vertex_pairs()
{
    return pattern_vertex_pairs;
}
const std::vector<bool> & Pattern::get_directions()
{
    return directions;
}
const std::vector<std::vector<timestamp_t>> & Pattern::get_timestamps()
{
    return timestamps;
}
const std::vector<std::vector<std::map<std::string, double>>> & Pattern::get_attributes()
{
    return attributes;
}

std::vector<std::vector<std::map<std::string, std::string>>> Pattern::get_decoded_attributes(std::map<std::string, AttributeType> edge_schema,
                                                                                              NominalEncoder * ne)
{

    std::vector<std::vector<std::map<std::string, std::string>>> pattern_attributes_decoded;
    for (int i = 0; i < attributes.size(); ++i)
    {
        std::vector<std::map<std::string, std::string>> edge_attributes;
        for (int j = 0; j < attributes.at(i).size(); ++j)
        {
            std::map<std::string, std::string> edge_instance_attributes;
            for (auto & kv : attributes.at(i).at(j))
            {
                if (edge_schema[kv.first] == AttributeType::NUMERIC)
                {
                    edge_instance_attributes[kv.first] = std::to_string(kv.second);
                }
                else
                {
                    edge_instance_attributes[kv.first] = ne->get_value(kv.second);
                }
            }
            edge_attributes.push_back(edge_instance_attributes);
        }
        pattern_attributes_decoded.push_back(edge_attributes);
    }
    return pattern_attributes_decoded;
}

const std::vector<double> & Pattern::get_scores()
{
    return scores;
}

void Pattern::clean_empty_instances()
{
    std::vector<std::vector<timestamp_t>> new_timestamps;
    std::vector<std::vector<std::map<std::string, double>>> new_attributes;

    for (int i = 0; i < attributes.size(); ++i)
    {
        std::vector<timestamp_t> new_timestamps_vector;
        std::vector<std::map<std::string, double>> new_attributes_vector;
        for (int j = 0; j < attributes[i].size(); ++j)
        {
            if (attributes[i][j].size() > 0)
            {
                new_timestamps_vector.push_back(timestamps[i][j]);
                new_attributes_vector.push_back(attributes[i][j]);
            }
        }
        new_timestamps.push_back(new_timestamps_vector);
        new_attributes.push_back(new_attributes_vector);
    }

    timestamps = new_timestamps;
    attributes = new_attributes;

}
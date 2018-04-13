#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <sstream>
#include "DynamicGraph.h"
#include "DynamicGraphExamples.h"
#include "RandomGenerator.h"
#include <ctime>

#include <random>
#include "PatternMiner.h"

#include "commonutilstemplated.h"


#include <experimental/filesystem>


namespace fs = std::experimental::filesystem;





template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


Vertex ver(int v_id, std::string label, NominalEncoder * ne)
{
    std::map<std::string, double> attributes;
    attributes["label"] = ne->get_encoding(label);
    return Vertex(v_id, attributes);
}

Edge edg(int e_id, int f_id, int t_id, std::string label, timestamp_t timestamp, NominalEncoder * ne)
{
    std::map<std::string, double> attributes;
    attributes["label"] = ne->get_encoding(label);
    return Edge(e_id, f_id, t_id, timestamp, attributes);
}

std::vector<Vertex> read_vertices(string directory, bool has_header, NominalEncoder * ne)
{
    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Vertex> vertices;

    for( std::string line; getline( input, line ); )
    {
        if (!has_header || line_number > 0)
        {
            std::vector<std::string> x = split(line, ',');
            vertices.push_back(ver(std::stoi(x[0]), x[1], ne));
        }
        line_number++;
    }
    return vertices;
}


std::vector<Edge> read_edges(string directory, bool has_header, NominalEncoder * ne)
{
    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Edge> edges;

    for( std::string line; getline( input, line ); )
    {
        if (!has_header || line_number > 0)
        {
            std::vector<std::string> x = split(line, ',');
            edges.push_back(edg(std::stoi(x[0]), std::stoi(x[1]), std::stoi(x[2]), x[3], std::stod(x[4]), ne));
        }
        line_number++;
    }
    return edges;
}


void output_base_edges_for_classifier(DynamicGraph * graph,
                                      std::vector<std::vector<int>> positive_event_edges,
                                      std::vector<std::vector<int>> positive_event_vertices,
                                      std::vector<timestamp_t > positive_event_times,
                                      std::vector<std::vector<int>> negative_event_edges,
                                      std::vector<std::vector<int>> negative_event_vertices,
                                      std::vector<timestamp_t > negative_event_times,
                                      std::string output_filename)
{

    ofstream myfile;
    myfile.open(output_filename);
    myfile << "id,label,timestamp,new_label,class" << std::endl;
    for (int i = 0; i < positive_event_vertices.size(); ++i)
    {
        for (int j = 0; j < positive_event_vertices[i].size(); ++j)
        {
            std::vector<Edge *> edges = graph->get_adjacency_list().at(positive_event_vertices[i][j]);
            for (auto & edge : edges)
            {
                // add only edges that are not the event edges
                // if there are no event edges, it is always ok
                if (positive_event_edges.size() == 0 || edge->get_original_edge_id() != positive_event_edges[0][j])
                {
//                    int time = (int) (positive_event_times[j] - edge->get_timestamp() + 1);
                    int time = (int) (positive_event_times[j] - edge->get_timestamp());
                    myfile << j << ","; // id
                    myfile << edge->get_attributes().at("label") << ","; // label
                    myfile << time << ","; // timestamp
                    myfile << "f" << edge->get_attributes().at("label") << "_" << time << ","; // new_label
                    myfile << "pos" << std::endl; // class
                }
            }
        }
    }

    for (int i = 0; i < negative_event_vertices.size(); ++i)
    {
        for (int j = 0; j < negative_event_vertices[i].size(); ++j)
        {
            std::vector<Edge *> edges = graph->get_adjacency_list().at(negative_event_vertices[i][j]);
            for (auto & edge : edges)
            {
                // add only edges that are not the event edges
                // if there are no event edges, it is always ok
                if (negative_event_edges.size() == 0 || edge->get_original_edge_id() != negative_event_edges[0][j])
                {
//                    int time = (int) (negative_event_times[j] - edge->get_timestamp() + 1);
                    int time = (int) (negative_event_times[j] - edge->get_timestamp());
                    myfile << (j + positive_event_times.size()) << ","; // id
                    myfile << edge->get_attributes().at("label") << ","; // label
                    myfile << time << ","; // timestamp
                    myfile << "f" << edge->get_attributes().at("label") << "_" << time << ","; // new_label
                    myfile << "neg" << std::endl; // class
                }
            }
        }
    }
    myfile.close();
}





void run_process(std::string vertex_filename, std::string edge_filename,
                 std::string positive_edge_label, std::string negative_edge_label,
                 std::string outputfile_prefix, timestamp_t time_unit_primary, timestamp_t time_unit_secondary,
                 bool undirected, int random_walks, double prob_restart, int max_pattern_edges, int evolution_epochs,
                 int evolution_subepochs, int evaluation_random_walks, bool use_vertex_attributes, int n_events,
                 int n_sample, int n_patterns, bool use_uniform_crossover, bool use_simple_init, bool limit_negative_population,
                 bool output_baseline_edges, unsigned seed,
                 bool verbose)
{

    fs::path p = fs::u8path(outputfile_prefix);
    if (!fs::is_directory(p) || !fs::exists(p)) { // Check if src folder exists
        fs::create_directories(p); // create src folder
    }

    clock_t begin_total = clock();
    clock_t begin_data_preparation = clock();

    RandomGenerator random_generator = RandomGenerator(seed);
//    RandomGenerator random_generator = RandomGenerator();
    NominalEncoder * ne = new NominalEncoder();

    if (verbose) print("LOADING THE GRAPH ... ");
    std::vector<Vertex> vertices = read_vertices(vertex_filename, true, ne);
    std::vector<Edge> edges = read_edges(edge_filename, true, ne);

    std::map<std::string, AttributeType> vertex_schema;
    std::map<std::string, AttributeType> edge_schema;
    vertex_schema["label"] = AttributeType::NOMINAL;
    edge_schema["label"] = AttributeType::NOMINAL;

    DynamicGraph graph = DynamicGraph(vertices, edges, vertex_schema, edge_schema, undirected);
    if (verbose) println("DONE");


    if (verbose) print("PREPARING EVENTS ... ");
    // find the events in the original graph:
    std::vector<int> all_positive_event_edges;
    std::vector<std::vector<int>> all_positive_event_vertices;
    std::vector<timestamp_t > all_positive_event_times;
    graph.identify_edge_events("label", (int) ne->get_encoding(positive_edge_label), n_events * 2, &random_generator,
                               all_positive_event_edges, all_positive_event_vertices, all_positive_event_times);
    std::vector<int> all_negative_event_edges;
    std::vector<std::vector<int>> all_negative_event_vertices;
    std::vector<timestamp_t > all_negative_event_times;
    graph.identify_edge_events("label", (int) ne->get_encoding(negative_edge_label), n_events * 2, &random_generator,
                               all_negative_event_edges, all_negative_event_vertices, all_negative_event_times);


    // divide the initial patterns into train and test sets:
    // train positive
    std::vector<std::vector<int>> positive_event_edges_train = {std::vector<int>()};
    std::vector<std::vector<int>> positive_event_vertices_train = {std::vector<int>(), std::vector<int>()};
    std::vector<timestamp_t > positive_event_times_train;
    // train negative
    std::vector<std::vector<int>> negative_event_edges_train = {std::vector<int>()};
    std::vector<std::vector<int>> negative_event_vertices_train = {std::vector<int>(), std::vector<int>()};
    std::vector<timestamp_t > negative_event_times_train;

    // test positive
    std::vector<std::vector<int>> positive_event_edges_test = {std::vector<int>()};
    std::vector<std::vector<int>> positive_event_vertices_test = {std::vector<int>(), std::vector<int>()};
    std::vector<timestamp_t > positive_event_times_test;
    // test negative
    std::vector<std::vector<int>> negative_event_edges_test = {std::vector<int>()};
    std::vector<std::vector<int>> negative_event_vertices_test = {std::vector<int>(), std::vector<int>()};
    std::vector<timestamp_t > negative_event_times_test;


    for (int i = 0; i < n_events * 2; ++i) {
        if (i < n_events)
        {
            // train set
            // positive
            positive_event_edges_train.at(0).push_back(all_positive_event_edges[i]);
            positive_event_vertices_train.at(0).push_back(all_positive_event_vertices[0][i]);
            positive_event_vertices_train.at(1).push_back(all_positive_event_vertices[1][i]);
            positive_event_times_train.push_back(all_positive_event_times[i]);
            // negative
            negative_event_edges_train.at(0).push_back(all_negative_event_edges[i]);
            negative_event_vertices_train.at(0).push_back(all_negative_event_vertices[0][i]);
            negative_event_vertices_train.at(1).push_back(all_negative_event_vertices[1][i]);
            negative_event_times_train.push_back(all_negative_event_times[i]);
        }
        else
        {
            // test set
            // positive
            positive_event_edges_test.at(0).push_back(all_positive_event_edges[i]);
            positive_event_vertices_test.at(0).push_back(all_positive_event_vertices[0][i]);
            positive_event_vertices_test.at(1).push_back(all_positive_event_vertices[1][i]);
            positive_event_times_test.push_back(all_positive_event_times[i]);
            // negative
            negative_event_edges_test.at(0).push_back(all_negative_event_edges[i]);
            negative_event_vertices_test.at(0).push_back(all_negative_event_vertices[0][i]);
            negative_event_vertices_test.at(1).push_back(all_negative_event_vertices[1][i]);
            negative_event_times_test.push_back(all_negative_event_times[i]);
        }
    }

    if (verbose) println("DONE");

    outputfile_prefix += "/";

    if (output_baseline_edges)
    {
        // print training edges
        output_base_edges_for_classifier(&graph,
                                         positive_event_edges_train, positive_event_vertices_train,
                                         positive_event_times_train, negative_event_edges_train,
                                         negative_event_vertices_train, negative_event_times_train,
                                         outputfile_prefix + "baseline_edges_data_train.csv");
        // print test edges
        output_base_edges_for_classifier(&graph,
                                         positive_event_edges_test, positive_event_vertices_test,
                                         positive_event_times_test, negative_event_edges_test,
                                         negative_event_vertices_test, negative_event_times_test,
                                         outputfile_prefix + "baseline_edges_data_test.csv");
    }

    std::vector<DynamicGraph> train_graph_instances = graph.create_subgraph_instances(positive_event_times_train,
                                                                                      negative_event_times_train,
                                                                                      time_unit_primary);

    std::vector<DynamicGraph> test_graph_instances = graph.create_subgraph_instances(positive_event_times_test,
                                                                                     negative_event_times_test,
                                                                                     time_unit_primary);

    clock_t end_data_preparation = clock();

    for (int pattern_number = 0; pattern_number < n_patterns; ++pattern_number)
    {
        if (verbose) println("------ SEARCHING FOR PATTERN No. ", pattern_number);

        clock_t begin_sampling = clock();

        // prepare a seed of events for pattern extraction
        // sample positive
        std::vector<std::vector<int>> positive_event_edges_sample = {std::vector<int>()};
        std::vector<std::vector<int>> positive_event_vertices_sample = {std::vector<int>(), std::vector<int>()};
        std::vector<timestamp_t > positive_event_times_sample;
        std::vector<int> selected_indices_positive = random_generator.generate_random_int_vector(n_sample,
                                                                                                 positive_event_times_train.size());
        for (int i = 0; i < n_sample; ++i)
        {
            positive_event_edges_sample.at(0).push_back(positive_event_edges_train[0][selected_indices_positive[i]]);
            positive_event_vertices_sample.at(0).push_back(positive_event_vertices_train[0][selected_indices_positive[i]]);
            positive_event_vertices_sample.at(1).push_back(positive_event_vertices_train[1][selected_indices_positive[i]]);
            positive_event_times_sample.push_back(positive_event_times_train[selected_indices_positive[i]]);
        }

        // sample negative
        std::vector<std::vector<int>> negative_event_edges_sample = {std::vector<int>()};
        std::vector<std::vector<int>> negative_event_vertices_sample = {std::vector<int>(), std::vector<int>()};
        std::vector<timestamp_t > negative_event_times_sample;
        std::vector<int> selected_indices_negative = random_generator.generate_random_int_vector(n_sample,
                                                                                                 negative_event_times_train.size());
        for (int i = 0; i < n_sample; ++i)
        {
            negative_event_edges_sample.at(0).push_back(negative_event_edges_train[0][selected_indices_negative[i]]);
            negative_event_vertices_sample.at(0).push_back(negative_event_vertices_train[0][selected_indices_negative[i]]);
            negative_event_vertices_sample.at(1).push_back(negative_event_vertices_train[1][selected_indices_negative[i]]);
            negative_event_times_sample.push_back(negative_event_times_train[selected_indices_negative[i]]);
        }
        clock_t end_sampling = clock();

        clock_t begin_pattern_mining = clock();
        std::vector<DynamicGraph> sample_graph_instances = graph.create_subgraph_instances(positive_event_times_sample,
                                                                                           negative_event_times_sample,
                                                                                           time_unit_primary);

        PatternMiner pattern_miner = PatternMiner(sample_graph_instances, positive_event_vertices_sample, positive_event_times_sample,
                                                  positive_event_edges_sample,
                                                  negative_event_vertices_sample, negative_event_times_sample,
                                                  negative_event_edges_sample,
                                                  use_vertex_attributes, time_unit_primary, time_unit_secondary,
                                                  random_walks, prob_restart, max_pattern_edges, evolution_epochs,
                                                  evolution_subepochs, &random_generator,
                                                  use_simple_init, use_uniform_crossover, limit_negative_population);

        std::vector<std::vector<double>> populations_fitness(max_pattern_edges * evolution_epochs, std::vector<double>());
        std::vector<std::vector<double>> negative_populations_fitness(max_pattern_edges * evolution_epochs * evolution_subepochs, std::vector<double>());


        Pattern pattern = pattern_miner.mine_pattern(populations_fitness, negative_populations_fitness, verbose);

        clock_t end_pattern_mining = clock();

        clock_t begin_evaluation = clock();
        if (verbose) print("ASSESSING PATTERN ... ");
        std::vector<std::vector<double>> evaluation_sample = pattern_miner.evaluate_pattern(&pattern, sample_graph_instances, positive_event_vertices_sample, positive_event_times_sample,
                                                                                            negative_event_vertices_sample, negative_event_times_sample, evaluation_random_walks);

        std::vector<std::vector<double>> evaluation_train = pattern_miner.evaluate_pattern(&pattern, train_graph_instances, positive_event_vertices_train, positive_event_times_train,
                                                                                           negative_event_vertices_train, negative_event_times_train, evaluation_random_walks);

        std::vector<std::vector<double>> evaluation_test = pattern_miner.evaluate_pattern(&pattern, test_graph_instances, positive_event_vertices_test, positive_event_times_test,
                                                                                          negative_event_vertices_test, negative_event_times_test, evaluation_random_walks);

        clock_t end_evaluation = clock();
        if (verbose) println("DONE");


        pattern.clean_empty_instances();

        double elapsed_secs_preparation = double(end_data_preparation - begin_data_preparation + end_sampling - begin_sampling) / CLOCKS_PER_SEC;
        double elapsed_secs_mining = double(end_pattern_mining - begin_pattern_mining) / CLOCKS_PER_SEC;
        double elapsed_secs_evaluation = double(end_evaluation - begin_evaluation) / CLOCKS_PER_SEC;

        std::string file_prefix = outputfile_prefix + "results_" + std::to_string(pattern_number);

        ofstream myfile;
        myfile.open(file_prefix + ".txt");
        // experiment parameters
        myfile << "DATA:" << outputfile_prefix << std::endl;
        myfile << "TOTAL_TIME_PREPARATION:" << std::to_string((int) elapsed_secs_preparation) << std::endl;
        myfile << "TOTAL_TIME_MINING:" << std::to_string((int) elapsed_secs_mining) << std::endl;
        myfile << "TOTAL_TIME_EVALUATION:" << std::to_string((int) elapsed_secs_evaluation) << std::endl;
        myfile << "TIME_UNIT_PRIMARY:" << std::to_string(time_unit_primary) << std::endl;
        myfile << "TIME_UNIT_SECONDARY:" << std::to_string(time_unit_secondary) << std::endl;
        myfile << "RANDOM_WALKS:" << std::to_string(random_walks) << std::endl;
        myfile << "PROB_RESTART:" << std::to_string(prob_restart) << std::endl;
        myfile << "MAX_PATTERN_EDGES:" << std::to_string(max_pattern_edges) << std::endl;
        myfile << "EVOLUTION_EPOCHS:" << std::to_string(evolution_epochs) << std::endl;
        myfile << "EVOLUTION_SUBEPOCHS:" << std::to_string(evolution_subepochs) << std::endl;
        myfile << "EVALUATION_RANDOM_WALKS:" << std::to_string(evaluation_random_walks) << std::endl;

        // pattern itself
        myfile << "PATTERN_EDGES:" << str(pattern.get_pattern_vertex_pairs()) << std::endl;
        myfile << "PATTERN_SCORES:" << str(pattern.get_scores()) << std::endl;
        myfile << "PATTERN_TIMESTAMPS:" << str(pattern.get_timestamps()) << std::endl;
        myfile << "PATTERN_ATTRIBUTES:" << str(pattern.get_decoded_attributes(graph.get_edge_schema(), ne)) << std::endl;
        myfile << "PATTERN_DIRECTIONS:" << str(pattern.get_directions()) << std::endl;
        myfile << "PATTERN_UNDIRECTED:" << str(graph.is_undirected()) << std::endl;

        // evaluation results
        myfile << "SAMPLE_EVALUATION_POSITIVE:" << str(evaluation_sample[0]) << std::endl;
        myfile << "SAMPLE_EVALUATION_NEGATIVE:" << str(evaluation_sample[1]) << std::endl;
        myfile << "TRAIN_EVALUATION_POSITIVE:" << str(evaluation_train[0]) << std::endl;
        myfile << "TRAIN_EVALUATION_NEGATIVE:" << str(evaluation_train[1]) << std::endl;
        myfile << "TEST_EVALUATION_POSITIVE:" << str(evaluation_test[0]) << std::endl;
        myfile << "TEST_EVALUATION_NEGATIVE:" << str(evaluation_test[1]) << std::endl;

        myfile.close();

        int number_of_quantiles = 5;
        ofstream myfile_positive;
        myfile_positive.open(file_prefix + "_positive_population.txt");
        for (int j = 0; j < populations_fitness.size(); ++j)
        {
            for (int i = 0; i < number_of_quantiles; ++i)
            {
                if (populations_fitness[j].size() > i)
                {
                    myfile_positive << populations_fitness[j][i];
                }
                if (i < number_of_quantiles - 1)
                {
                    myfile_positive << ",";
                }
            }
            myfile_positive << std::endl;
        }
        myfile_positive.close();


        ofstream myfile_negative;
        myfile_negative.open(file_prefix + "_negative_population.txt");
        for (int j = 0; j < negative_populations_fitness.size(); ++j)
        {
            for (int i = 0; i < number_of_quantiles; ++i)
            {
                if (negative_populations_fitness[j].size() > i)
                {
                    myfile_negative << negative_populations_fitness[j][i];
                }
                if (i < number_of_quantiles - 1)
                {
                    myfile_negative << ",";
                }
            }
            myfile_negative << std::endl;
        }
        myfile_negative.close();

    }

    delete ne;

    clock_t end_total = clock();

    double elapsed_secs_total = double(end_total - begin_total) / CLOCKS_PER_SEC;
    if (verbose) println("EXPERIMENT_TOTAL_TIME: ", ((int) elapsed_secs_total));
}





void run_dblp_experiment(int experiment)
{
    std::vector<std::string> outputfile_prefix_params = {"results/dblp1", "results/dblp2", "results/dblp3", "results/dblp4"};
    std::vector<std::string> positive_edge_label_params = {"icml", "kdd", "nips", "nips"};
    std::vector<std::string> negative_edge_label_params = {"kdd", "icml", "kdd", "icml"};

    std::string vertex_filename = "data/dblp_vertices.csv";
    std::string edge_filename = "data/dblp_edges.csv";

    std::string outputfile_prefix = outputfile_prefix_params[experiment];

    std::string positive_edge_label = positive_edge_label_params[experiment];
    std::string negative_edge_label = negative_edge_label_params[experiment];

    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 0.5;
    bool undirected = true;
    int random_walks = 1000;
    double prob_restart = 0.3;

    int max_pattern_edges = 5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_vertex_attributes = false;
    unsigned seed = 1;
    // how many events from positive (resp. negative) set into training (resp. testing) set
    int n_events = 100;
    // how many positive (resp. negative) events from training set are used for extracting a pattern
    int n_sample = 30;
    // how many patterns to extract, i.e. how many samplings to use
    int n_patterns = 20;

    bool use_uniform_crossover = true;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    bool output_baseline_edges = true;
    bool verbose = true;


    run_process(vertex_filename, edge_filename, positive_edge_label, negative_edge_label,
                outputfile_prefix, time_unit_primary, time_unit_secondary,
                undirected, random_walks, prob_restart, max_pattern_edges, evolution_epochs,
                evolution_subepochs, evaluation_random_walks, use_vertex_attributes, n_events,
                n_sample, n_patterns, use_uniform_crossover, use_simple_init, limit_negative_population,
                output_baseline_edges, seed,
                verbose);
}





void run_enron_experiment()
{

    std::string vertex_filename = "data/enron_vertices.csv";
    std::string edge_filename = "data/enron_edges.csv";

    std::string outputfile_prefix = "results/enron";

    std::string positive_edge_label = "2";
    std::string negative_edge_label = "5";

    timestamp_t time_unit_primary = 10 * 48 * 3600.0;
    timestamp_t time_unit_secondary = 2 * 48 * 3600.0;
    bool undirected = false;
    int random_walks = 3000;
    double prob_restart = 0.3;

    int max_pattern_edges = 5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_vertex_attributes = false;
    unsigned seed = 1;
    // how many events from positive (resp. negative) set into training (resp. testing) set
    int n_events = 20;
    // how many positive (resp. negative) events from training set are used for extracting a pattern
    int n_sample = 10;
    // how many patterns to extract, i.e. how many samplings to use
    int n_patterns = 10;

    bool use_uniform_crossover = true;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    bool output_baseline_edges = true;
    bool verbose = true;


    run_process(vertex_filename, edge_filename, positive_edge_label, negative_edge_label,
                outputfile_prefix, time_unit_primary, time_unit_secondary,
                undirected, random_walks, prob_restart, max_pattern_edges, evolution_epochs,
                evolution_subepochs, evaluation_random_walks, use_vertex_attributes, n_events,
                n_sample, n_patterns, use_uniform_crossover, use_simple_init, limit_negative_population,
                output_baseline_edges, seed,
                verbose);

}



void print_usage()
{
    std::cerr << "EWALDIS takes the following arguments only: dblp1, dblp2, dblp3, dblp4, enron" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        print_usage();
    }

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "dblp1") run_dblp_experiment(0);
        else if (arg == "dblp2") run_dblp_experiment(1);
        else if (arg == "dblp3") run_dblp_experiment(2);
        else if (arg == "dblp4") run_dblp_experiment(3);
        else if (arg == "enron") run_enron_experiment();
        else {
            std::cerr << "Cannot recognize argument '" << arg << "'." << std::endl;
            print_usage();
        }
    }

    return 0;
}






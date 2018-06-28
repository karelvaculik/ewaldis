#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <sstream>
#include "ewaldis/DynamicGraph.h"
#include "ewaldis/DynamicGraphExamples.h"
#include "ewaldis/RandomGenerator.h"
#include <ctime>

#include <random>
#include "PatternMiner.h"

#include "commonutilstemplated.h"


#include <experimental/filesystem>


namespace fs = std::experimental::filesystem;


static void show_usage(std::string name) {
    std::cerr << "EWALDIS, version 1.0.0"
              << std::endl
              << "Author: Karel Vaculik (vaculik.dev@gmail.com)"
              << std::endl
              << std::endl
              << "Usage: " << name << " <option(s)>"
              << std::endl
              << std::endl
              << "NOTE: N stands for an integer, R for a real number, F a filename string"
              << std::endl
              << std::endl
              << "\t-h, --help\t\t\tShows this help message"
              << std::endl
              << "\t-v, --verbose\t\tVerbose output"
              << std::endl
              << std::endl
              << "MANDATORY ARGUMENTS:"
              << std::endl
              << "\t--vertices F\t\tSpecifies the input file with vertices"
              << std::endl
              << "\t--edges F\t\t\tSpecifies the input file with edges"
              << std::endl
              << "\t--train_pos F\t\tSpecifies the training input file with positive events"
              << std::endl
              << "\t--train_neg F\t\tSpecifies the training input file with negative events"
              << std::endl
              << "\t-o OUTPUT_DIR\t\tSpecifies the output directory"
              << std::endl
              << "\t-e {vertices,edges}\tSpecifies the type of events: either vertices or edges"
              << std::endl
              << "\t-a {nominal,numerical}\tSpecifies the type of edge attributes: either nominal or numerical"
              << std::endl
              << std::endl
              << "OPTIONAL ARGUMENTS:"
              << std::endl
              << "\t-p N\t\t\t\tSpecifies the number of patterns to be mined; 1 by default"
              << std::endl
              << "\t-n N\t\t\t\tSpecifies the size of sample from positive and negative events used for pattern mining; if not specified, all events are used"
              << std::endl
              << "\t--test_pos F\t\tSpecifies the test input file with positive events"
              << std::endl
              << "\t--test_neg F\t\tSpecifies the test input file with negative events"
              << std::endl
              << "\t-u\t\t\t\t\tTreat the input graph as undirected; if not specified, it is considered directed"
              << std::endl
              << "\t-m N\t\t\t\tSpecifies the maximum number of edges per pattern; 20 by default"
              << std::endl
              << "\t-w N\t\t\t\tSpecifies the number of random walks; 1000 by default"
              << std::endl
              << "\t-r R\t\t\t\tSpecifies the probability of random walks restart; decimal value from [0.0, 1.0]; 0.3 by default"
              << std::endl
              << "\t-pu R\t\t\t\tSpecifies the value of primary time unit; 1.0 by default"
              << std::endl
              << "\t-su R\t\t\t\tSpecifies the value of secondary time unit; 1.0 by default"
              << std::endl
              << "\t-ep N\t\t\t\tSpecifies the number of epochs in the genetic algorithm; 25 by default"
              << std::endl
              << "\t-se N\t\t\t\tSpecifies the number of subepochs in the genetic algorithm; 25 by default"
              << std::endl
              << "\t-ew N\t\t\t\tSpecifies the number of random walks used for evaluation; 10 by default"
              << std::endl
              << "\t--unicross\t\t\tUse uniform crossover instead of single-point crossover"
              << std::endl
              << "\t--baseline F\t\tOutput edges for baseline classifier"
              << std::endl
              << "\t-s N\t\t\t\tSpecifies seed for random generator"
              << std::endl;
}


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


Vertex ver(int v_id, std::string label, NominalEncoder *ne) {
    std::map<std::string, double> attributes;
    attributes["label"] = ne->get_encoding(label);
    return Vertex(v_id, attributes);
}

Edge edg(int e_id, int f_id, int t_id, std::string label, timestamp_t timestamp, NominalEncoder *ne) {
    std::map<std::string, double> attributes;
    attributes["label"] = ne->get_encoding(label);
    return Edge(e_id, f_id, t_id, timestamp, attributes);
}

std::vector<Vertex> read_vertices(string directory, bool has_header, NominalEncoder *ne) {
    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Vertex> vertices;

    for (std::string line; getline(input, line);) {
        if (!has_header || line_number > 0) {
            std::vector<std::string> x = split(line, ',');
            vertices.push_back(ver(std::stoi(x[0]), x[1], ne));
        }
        line_number++;
    }
    return vertices;
}


std::vector<Edge> read_edges(string directory, bool has_header, NominalEncoder *ne) {
    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Edge> edges;

    for (std::string line; getline(input, line);) {
        if (!has_header || line_number > 0) {
            std::vector<std::string> x = split(line, ',');
            edges.push_back(edg(std::stoi(x[0]), std::stoi(x[1]), std::stoi(x[2]), x[3], std::stod(x[4]), ne));
        }
        line_number++;
    }
    return edges;
}

void read_vertex_events(std::string directory, bool has_header,
                        std::vector<std::vector<int>> &vertex_ids,
                        std::vector<timestamp_t> &vertex_timestamps) {
    vertex_ids.push_back(std::vector<int>());

    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Edge> edges;

    for (std::string line; getline(input, line);) {
        if (!has_header || line_number > 0) {
            std::vector<std::string> x = split(line, ',');
            vertex_ids[0].push_back(std::stoi(x[0]));
            vertex_timestamps.push_back(std::stod(x[1]));
        }
        line_number++;
    }
}

void read_edge_events(std::string directory, bool has_header,
                      std::vector<std::vector<int>> &edge_ids,
                      std::vector<std::vector<int>> &vertex_ids,
                      std::vector<timestamp_t> &vertex_timestamps) {
    edge_ids.push_back(std::vector<int>());
    vertex_ids.push_back(std::vector<int>());
    vertex_ids.push_back(std::vector<int>());

    std::ifstream input(directory);
    int line_number = 0;
    std::vector<Edge> edges;

    for (std::string line; getline(input, line);) {
        if (!has_header || line_number > 0) {
            std::vector<std::string> x = split(line, ',');
            edge_ids[0].push_back(std::stoi(x[0]));
            vertex_ids[0].push_back(std::stoi(x[1]));
            vertex_ids[1].push_back(std::stoi(x[2]));
            vertex_timestamps.push_back(std::stod(x[3]));
        }
        line_number++;
    }
}


void output_base_edges_for_classifier(DynamicGraph *graph,
                                      std::vector<std::vector<int>> positive_event_edges,
                                      std::vector<std::vector<int>> positive_event_vertices,
                                      std::vector<timestamp_t> positive_event_times,
                                      std::vector<std::vector<int>> negative_event_edges,
                                      std::vector<std::vector<int>> negative_event_vertices,
                                      std::vector<timestamp_t> negative_event_times,
                                      std::string output_filename) {

    ofstream myfile;
    myfile.open(output_filename);
    myfile << "id,label,timestamp,new_label,class" << std::endl;
    for (int i = 0; i < positive_event_vertices.size(); ++i) {
        for (int j = 0; j < positive_event_vertices[i].size(); ++j) {
            std::vector<Edge *> edges = graph->get_adjacency_list().at(positive_event_vertices[i][j]);
            for (auto &edge : edges) {
                // add only edges that are not the event edges
                // if there are no event edges, it is always ok
                if (positive_event_edges.size() == 0 || edge->get_original_edge_id() != positive_event_edges[0][j]) {
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

    for (int i = 0; i < negative_event_vertices.size(); ++i) {
        for (int j = 0; j < negative_event_vertices[i].size(); ++j) {
            std::vector<Edge *> edges = graph->get_adjacency_list().at(negative_event_vertices[i][j]);
            for (auto &edge : edges) {
                // add only edges that are not the event edges
                // if there are no event edges, it is always ok
                if (negative_event_edges.size() == 0 || edge->get_original_edge_id() != negative_event_edges[0][j]) {
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


void run_process_new(std::string vertex_filename, std::string edge_filename,
                     std::string train_pos_filename, std::string train_neg_filename,
                     std::string output_dir, bool vertex_events, std::string edge_attributes_type,
                     int n_patterns, bool use_sample, int n_sample,
                     std::string test_pos_filename, std::string test_neg_filename,
                     bool undirected, int max_pattern_edges, int random_walks, double prob_restart,
                     timestamp_t time_unit_primary, timestamp_t time_unit_secondary,
                     int evolution_epochs, int evolution_subepochs, int evaluation_random_walks,
                     bool use_uniform_crossover, bool output_baseline_edges, bool use_seed, unsigned seed,
                     bool use_vertex_attributes, bool use_simple_init, bool limit_negative_population,
                     bool verbose) {

    clock_t begin_total = clock();
    clock_t begin_data_preparation = clock();

    RandomGenerator random_generator;
    if (use_seed) {
        random_generator = RandomGenerator(seed);
    } else {
        random_generator = RandomGenerator();
    }
    NominalEncoder *ne = new NominalEncoder();

    if (verbose) print("LOADING THE GRAPH ... ");
    std::vector<Vertex> vertices = read_vertices(vertex_filename, true, ne);
    std::vector<Edge> edges = read_edges(edge_filename, true, ne);

    std::map<std::string, AttributeType> vertex_schema;
    std::map<std::string, AttributeType> edge_schema;
    vertex_schema["label"] = AttributeType::NOMINAL;
    if (edge_attributes_type == "nominal") {
        edge_schema["label"] = AttributeType::NOMINAL;
    } else {
        edge_schema["label"] = AttributeType::NUMERIC;
    }


    DynamicGraph graph = DynamicGraph(vertices, edges, vertex_schema, edge_schema, undirected);
    if (verbose) println("DONE");


    if (verbose) print("PREPARING EVENTS ... ");

    // train positive
    std::vector<std::vector<int>> positive_event_edges_train;
    std::vector<std::vector<int>> positive_event_vertices_train;
    std::vector<timestamp_t> positive_event_times_train;
    // train negative
    std::vector<std::vector<int>> negative_event_edges_train;
    std::vector<std::vector<int>> negative_event_vertices_train;
    std::vector<timestamp_t> negative_event_times_train;

    // test positive
    std::vector<std::vector<int>> positive_event_edges_test;
    std::vector<std::vector<int>> positive_event_vertices_test;
    std::vector<timestamp_t> positive_event_times_test;
    // test negative
    std::vector<std::vector<int>> negative_event_edges_test;
    std::vector<std::vector<int>> negative_event_vertices_test;
    std::vector<timestamp_t> negative_event_times_test;


    if (vertex_events) {
        read_vertex_events(train_pos_filename, true, positive_event_vertices_train, positive_event_times_train);
        read_vertex_events(train_neg_filename, true, negative_event_vertices_train, negative_event_times_train);

        if (test_pos_filename != "") {
            read_vertex_events(test_pos_filename, true, positive_event_vertices_test, positive_event_times_test);
            read_vertex_events(test_neg_filename, true, negative_event_vertices_test, negative_event_times_test);
        }
    } else {
        read_edge_events(train_pos_filename, true, positive_event_edges_train, positive_event_vertices_train,
                         positive_event_times_train);
        read_edge_events(train_neg_filename, true, negative_event_edges_train, negative_event_vertices_train,
                         negative_event_times_train);

        if (test_pos_filename != "") {
            read_edge_events(test_pos_filename, true, positive_event_edges_test, positive_event_vertices_test,
                             positive_event_times_test);
            read_edge_events(test_neg_filename, true, negative_event_edges_test, negative_event_vertices_test,
                             negative_event_times_test);
        }
    }

    graph.get_edges().at()

    return ;

    // create directory if it does not exist
    fs::path p = fs::u8path(output_dir);
    if (!fs::is_directory(p) || !fs::exists(p)) { // Check if src folder exists
        fs::create_directories(p); // create src folder
    }

    if (verbose) println("DONE");

    if (output_baseline_edges) {
        // print training edges
        std::cout << "PRINTING BASELINE";
        output_base_edges_for_classifier(&graph,
                                         positive_event_edges_train, positive_event_vertices_train,
                                         positive_event_times_train, negative_event_edges_train,
                                         negative_event_vertices_train, negative_event_times_train,
                                         output_dir + "/baseline_edges_data_train.csv");
        // print test edges
        output_base_edges_for_classifier(&graph,
                                         positive_event_edges_test, positive_event_vertices_test,
                                         positive_event_times_test, negative_event_edges_test,
                                         negative_event_vertices_test, negative_event_times_test,
                                         output_dir + "/baseline_edges_data_test.csv");
    }

    return ;

    std::vector<DynamicGraph> train_graph_instances = graph.create_subgraph_instances(positive_event_times_train,
                                                                                      negative_event_times_train,
                                                                                      time_unit_primary);

    std::vector<DynamicGraph> test_graph_instances = graph.create_subgraph_instances(positive_event_times_test,
                                                                                     negative_event_times_test,
                                                                                     time_unit_primary);

    clock_t end_data_preparation = clock();

    for (int pattern_number = 0; pattern_number < n_patterns; ++pattern_number) {
        if (verbose) println("------ SEARCHING FOR PATTERN No. ", pattern_number);

        clock_t begin_sampling = clock();

        // prepare a seed of events for pattern extraction
        // sample positive
        std::vector<std::vector<int>> positive_event_edges_sample = {std::vector<int>()};
        std::vector<std::vector<int>> positive_event_vertices_sample = {std::vector<int>(), std::vector<int>()};
        std::vector<timestamp_t> positive_event_times_sample;
        std::vector<int> selected_indices_positive = random_generator.generate_random_int_vector(n_sample,
                                                                                                 positive_event_times_train.size());
        for (int i = 0; i < n_sample; ++i) {
            positive_event_edges_sample.at(0).push_back(positive_event_edges_train[0][selected_indices_positive[i]]);
            positive_event_vertices_sample.at(0).push_back(
                    positive_event_vertices_train[0][selected_indices_positive[i]]);
            positive_event_vertices_sample.at(1).push_back(
                    positive_event_vertices_train[1][selected_indices_positive[i]]);
            positive_event_times_sample.push_back(positive_event_times_train[selected_indices_positive[i]]);
        }

        // sample negative
        std::vector<std::vector<int>> negative_event_edges_sample = {std::vector<int>()};
        std::vector<std::vector<int>> negative_event_vertices_sample = {std::vector<int>(), std::vector<int>()};
        std::vector<timestamp_t> negative_event_times_sample;
        std::vector<int> selected_indices_negative = random_generator.generate_random_int_vector(n_sample,
                                                                                                 negative_event_times_train.size());
        for (int i = 0; i < n_sample; ++i) {
            negative_event_edges_sample.at(0).push_back(negative_event_edges_train[0][selected_indices_negative[i]]);
            negative_event_vertices_sample.at(0).push_back(
                    negative_event_vertices_train[0][selected_indices_negative[i]]);
            negative_event_vertices_sample.at(1).push_back(
                    negative_event_vertices_train[1][selected_indices_negative[i]]);
            negative_event_times_sample.push_back(negative_event_times_train[selected_indices_negative[i]]);
        }
        clock_t end_sampling = clock();

        clock_t begin_pattern_mining = clock();
        std::vector<DynamicGraph> sample_graph_instances = graph.create_subgraph_instances(positive_event_times_sample,
                                                                                           negative_event_times_sample,
                                                                                           time_unit_primary);

        PatternMiner pattern_miner = PatternMiner(sample_graph_instances, positive_event_vertices_sample,
                                                  positive_event_times_sample,
//        PatternMiner pattern_miner = PatternMiner(&graph, positive_event_vertices_sample, positive_event_times_sample,
                                                  positive_event_edges_sample,
                                                  negative_event_vertices_sample, negative_event_times_sample,
                                                  negative_event_edges_sample,
                                                  use_vertex_attributes, time_unit_primary, time_unit_secondary,
                                                  random_walks, prob_restart, max_pattern_edges, evolution_epochs,
                                                  evolution_subepochs, &random_generator,
                                                  use_simple_init, use_uniform_crossover, limit_negative_population);

        std::vector<std::vector<double>> populations_fitness(max_pattern_edges * evolution_epochs,
                                                             std::vector<double>());
        std::vector<std::vector<double>> negative_populations_fitness(
                max_pattern_edges * evolution_epochs * evolution_subepochs, std::vector<double>());


        Pattern pattern = pattern_miner.mine_pattern(populations_fitness, negative_populations_fitness, verbose);

        clock_t end_pattern_mining = clock();

        clock_t begin_evaluation = clock();
        if (verbose) print("ASSESSING PATTERN ... ");
        std::vector<std::vector<double>> evaluation_sample = pattern_miner.evaluate_pattern(&pattern,
                                                                                            sample_graph_instances,
                                                                                            positive_event_vertices_sample,
                                                                                            positive_event_times_sample,
                                                                                            negative_event_vertices_sample,
                                                                                            negative_event_times_sample,
                                                                                            evaluation_random_walks);

        std::vector<std::vector<double>> evaluation_train = pattern_miner.evaluate_pattern(&pattern,
                                                                                           train_graph_instances,
                                                                                           positive_event_vertices_train,
                                                                                           positive_event_times_train,
                                                                                           negative_event_vertices_train,
                                                                                           negative_event_times_train,
                                                                                           evaluation_random_walks);

        std::vector<std::vector<double>> evaluation_test = pattern_miner.evaluate_pattern(&pattern,
                                                                                          test_graph_instances,
                                                                                          positive_event_vertices_test,
                                                                                          positive_event_times_test,
                                                                                          negative_event_vertices_test,
                                                                                          negative_event_times_test,
                                                                                          evaluation_random_walks);
//        std::vector<std::vector<double>> evaluation_sample = pattern_miner.evaluate_pattern(&pattern, &graph, positive_event_vertices_sample, positive_event_times_sample,
//                                                                                            negative_event_vertices_sample, negative_event_times_sample, evaluation_random_walks);
//
//        std::vector<std::vector<double>> evaluation_train = pattern_miner.evaluate_pattern(&pattern, &graph, positive_event_vertices_train, positive_event_times_train,
//                                                                                           negative_event_vertices_train, negative_event_times_train, evaluation_random_walks);
//
//        std::vector<std::vector<double>> evaluation_test = pattern_miner.evaluate_pattern(&pattern, &graph, positive_event_vertices_test, positive_event_times_test,
//                                                                                          negative_event_vertices_test, negative_event_times_test, evaluation_random_walks);

        clock_t end_evaluation = clock();
        if (verbose) println("DONE");


        pattern.clean_empty_instances();

        double elapsed_secs_preparation =
                double(end_data_preparation - begin_data_preparation + end_sampling - begin_sampling) / CLOCKS_PER_SEC;
        double elapsed_secs_mining = double(end_pattern_mining - begin_pattern_mining) / CLOCKS_PER_SEC;
        double elapsed_secs_evaluation = double(end_evaluation - begin_evaluation) / CLOCKS_PER_SEC;

        std::string file_prefix = output_dir + "/results_" + std::to_string(pattern_number);

        ofstream myfile;
        myfile.open(file_prefix + ".txt");
        // experiment parameters
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
        myfile << "PATTERN_ATTRIBUTES:" << str(pattern.get_decoded_attributes(graph.get_edge_schema(), ne))
               << std::endl;
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
        for (int j = 0; j < populations_fitness.size(); ++j) {
            for (int i = 0; i < number_of_quantiles; ++i) {
                if (populations_fitness[j].size() > i) {
                    myfile_positive << populations_fitness[j][i];
                }
                if (i < number_of_quantiles - 1) {
                    myfile_positive << ",";
                }
            }
            myfile_positive << std::endl;
        }
        myfile_positive.close();


        ofstream myfile_negative;
        myfile_negative.open(file_prefix + "_negative_population.txt");
        for (int j = 0; j < negative_populations_fitness.size(); ++j) {
            for (int i = 0; i < number_of_quantiles; ++i) {
                if (negative_populations_fitness[j].size() > i) {
                    myfile_negative << negative_populations_fitness[j][i];
                }
                if (i < number_of_quantiles - 1) {
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


void run_dblp_1() {
    bool verbose = true;

    // mandatory arguments
    std::string vertex_filename = "data/dblp1/dblp_vertices.csv";
    std::string edge_filename = "data/dblp1/dblp_edges.csv";
    std::string train_pos_filename = "data/dblp1/dblp_train_positive_events.csv";
    std::string train_neg_filename = "data/dblp1/dblp_train_negative_events.csv";
    std::string output_dir = "results/dblp1_new3";
    bool vertex_events = false;
    std::string edge_attributes_type = "nominal";

    // optional arguments
    int n_patterns = 20;
    bool use_sample = true;
    int n_sample = 30;
    std::string test_pos_filename = "data/dblp1/dblp_test_positive_events.csv";
    std::string test_neg_filename = "data/dblp1/dblp_test_negative_events.csv";
    bool undirected = true;
    int max_pattern_edges = 5;
    int random_walks = 1000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 0.5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = true;
    bool output_baseline_edges = true;
    bool use_seed = true;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);
}


void run_dblp_2() {
    bool verbose = true;

    // mandatory arguments
    std::string vertex_filename = "data/dblp2/dblp_vertices.csv";
    std::string edge_filename = "data/dblp2/dblp_edges.csv";
    std::string train_pos_filename = "data/dblp2/dblp_train_positive_events.csv";
    std::string train_neg_filename = "data/dblp2/dblp_train_negative_events.csv";
    std::string output_dir = "results/dblp2_new3";
    bool vertex_events = false;
    std::string edge_attributes_type = "nominal";

    // optional arguments
    int n_patterns = 20;
    bool use_sample = true;
    int n_sample = 30;
    std::string test_pos_filename = "data/dblp2/dblp_test_positive_events.csv";
    std::string test_neg_filename = "data/dblp2/dblp_test_negative_events.csv";
    bool undirected = true;
    int max_pattern_edges = 5;
    int random_walks = 1000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 0.5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = true;
    bool output_baseline_edges = true;
    bool use_seed = true;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);
}


void run_dblp_3() {
    bool verbose = true;

    // mandatory arguments
    std::string vertex_filename = "data/dblp3/dblp_vertices.csv";
    std::string edge_filename = "data/dblp3/dblp_edges.csv";
    std::string train_pos_filename = "data/dblp3/dblp_train_positive_events.csv";
    std::string train_neg_filename = "data/dblp3/dblp_train_negative_events.csv";
    std::string output_dir = "results/dblp3_new2";
    bool vertex_events = false;
    std::string edge_attributes_type = "nominal";

    // optional arguments
    int n_patterns = 20;
    bool use_sample = true;
    int n_sample = 30;
    std::string test_pos_filename = "data/dblp3/dblp_test_positive_events.csv";
    std::string test_neg_filename = "data/dblp3/dblp_test_negative_events.csv";
    bool undirected = true;
    int max_pattern_edges = 5;
    int random_walks = 1000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 0.5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = true;
    bool output_baseline_edges = true;
    bool use_seed = true;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);
}


void run_dblp_4() {
    bool verbose = true;

    // mandatory arguments
    std::string vertex_filename = "data/dblp4/dblp_vertices.csv";
    std::string edge_filename = "data/dblp4/dblp_edges.csv";
    std::string train_pos_filename = "data/dblp4/dblp_train_positive_events.csv";
    std::string train_neg_filename = "data/dblp4/dblp_train_negative_events.csv";
    std::string output_dir = "results/dblp4_new2";
    bool vertex_events = false;
    std::string edge_attributes_type = "nominal";

    // optional arguments
    int n_patterns = 20;
    bool use_sample = true;
    int n_sample = 30;
    std::string test_pos_filename = "data/dblp4/dblp_test_positive_events.csv";
    std::string test_neg_filename = "data/dblp4/dblp_test_negative_events.csv";
    bool undirected = true;
    int max_pattern_edges = 5;
    int random_walks = 1000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 0.5;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = true;
    bool output_baseline_edges = true;
    bool use_seed = true;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);
}


void run_enron() {
    bool verbose = true;

    // mandatory arguments
    std::string vertex_filename = "data/enron/enron_vertices.csv";
    std::string edge_filename = "data/enron/enron_edges.csv";
    std::string train_pos_filename = "data/enron/enron_train_positive_events.csv";
    std::string train_neg_filename = "data/enron/enron_train_negative_events.csv";
    std::string output_dir = "results/enron_new2";
    bool vertex_events = false;
    std::string edge_attributes_type = "nominal";

    // optional arguments
    int n_patterns = 10;
    bool use_sample = true;
    int n_sample = 10;
    std::string test_pos_filename = "data/enron/enron_test_positive_events.csv";
    std::string test_neg_filename = "data/enron/enron_test_negative_events.csv";
    bool undirected = false;
    int max_pattern_edges = 5;
    int random_walks = 3000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 10 * 48 * 3600.0;
    timestamp_t time_unit_secondary = 2 * 48 * 3600.0;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = true;
    bool output_baseline_edges = true;
    bool use_seed = true;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);
}


int main(int argc, char *argv[]) {
//    run_dblp_1();
//    run_dblp_2();
//    run_dblp_3();
//    run_dblp_4();
//    run_enron();


//    return 0;

    bool verbose = false;

    // mandatory arguments
    std::string vertex_filename = "";
    std::string edge_filename = "";
    std::string train_pos_filename = "";
    std::string train_neg_filename = "";
    std::string output_dir = "";
    bool vertex_events;
    std::string str_vertex_events = "";
    std::string edge_attributes_type = "";

    // optional arguments
    int n_patterns = 1;
    bool use_sample = false;
    int n_sample = 0;
    std::string test_pos_filename = "";
    std::string test_neg_filename = "";
    bool undirected = false;
    int max_pattern_edges = 20;
    int random_walks = 1000;
    double prob_restart = 0.3;
    timestamp_t time_unit_primary = 1.0;
    timestamp_t time_unit_secondary = 1.0;
    int evolution_epochs = 25;
    int evolution_subepochs = 25;
    int evaluation_random_walks = 10;
    bool use_uniform_crossover = false;
    bool output_baseline_edges = false;
    bool use_seed = false;
    unsigned seed = 1;


    // these cannot be changed at this moment
    bool use_vertex_attributes = false;
    bool use_simple_init = false;
    bool limit_negative_population = true;

    // parse the arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-v") || (arg == "--verbose")) {
            verbose = true;
        }
            // STRINGS
        else if (arg == "--vertices" || arg == "--edges" || arg == "--train_pos" || arg == "--train_neg" ||
                 arg == "-o" || arg == "--test_pos" || arg == "--test_neg" || arg == "--baseline") {
            if (i + 1 >= argc) {
                std::cerr << arg << " option requires one argument." << std::endl;
                return 1;
            } else if (arg == "--vertices") {
                vertex_filename = argv[++i];
            } else if (arg == "--edges") {
                edge_filename = argv[++i];
            } else if (arg == "--train_pos") {
                train_pos_filename = argv[++i];
            } else if (arg == "--train_neg") {
                train_neg_filename = argv[++i];
            } else if (arg == "-o") {
                output_dir = argv[++i];
            } else if (arg == "--test_pos") {
                test_pos_filename = argv[++i];
            } else if (arg == "--test_neg") {
                test_neg_filename = argv[++i];
            }
        }
            // INTEGERS
        else if (arg == "-p" || arg == "-n" || arg == "-m" || arg == "-w" || arg == "-ep" || arg == "-se" ||
                 arg == "-ew" || arg == "-s") {
            if (i + 1 >= argc) {
                std::cerr << arg << " option requires one argument." << std::endl;
                return 1;
            }
            std::string str_number = argv[++i];
            try {
                int number = std::stoi(str_number);
                if (arg == "-p") {
                    n_patterns = number;
                } else if (arg == "-n") {
                    n_sample = number;
                    if (n_sample > 0) use_sample = true;
                } else if (arg == "-m") {
                    max_pattern_edges = number;
                } else if (arg == "-w") {
                    random_walks = number;
                } else if (arg == "-ep") {
                    evolution_epochs = number;
                } else if (arg == "-se") {
                    evolution_subepochs = number;
                } else if (arg == "-ew") {
                    evaluation_random_walks = number;
                } else if (arg == "-s") {
                    seed = number;
                    use_seed = true;
                }

            }
            catch (std::exception &) {
                std::cerr << arg << " option requires an integer argument." << std::endl;
            }

        }
            // REALS
        else if (arg == "-r" || arg == "-pu" || arg == "-su") {
            if (i + 1 >= argc) {
                std::cerr << arg << " option requires one argument." << std::endl;
                return 1;
            }
            std::string str_number = argv[++i];
            try {
                double number = std::stod(str_number);
                if (arg == "-r") {
                    if (number < 0 || number > 1) {
                        std::cerr << "-r option requires a numeric argument from interval [0, 1]" << std::endl;
                        return 1;
                    }
                    prob_restart = number;
                } else if (arg == "-pu") {
                    if (number <= 0) {
                        std::cerr << "-pu option requires a positive numeric argument" << std::endl;
                        return 1;
                    }
                    time_unit_primary = number;
                } else if (arg == "-su") {
                    if (number <= 0) {
                        std::cerr << "-su option requires a positive numeric argument" << std::endl;
                        return 1;
                    }
                    time_unit_secondary = number;
                }
            }
            catch (std::exception &) {
                std::cerr << arg << " option requires a numeric argument." << std::endl;
                return 1;
            }

        }
            // OTHERS
        else if (arg == "-e") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                str_vertex_events = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
                if (str_vertex_events != "vertices" && str_vertex_events != "edges") {
                    std::cerr << "-e option requires either vertices or edges." << std::endl;
                    return 1;
                } else {
                    vertex_events = str_vertex_events == "vertices";
                }
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-e option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-a") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                edge_attributes_type = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
                if (edge_attributes_type != "nominal" && edge_attributes_type != "numerical") {
                    std::cerr << "-a option requires either nominal or numerical." << std::endl;
                    return 1;
                }
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-a option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-u") {
            undirected = true;
        } else if (arg == "--unicross") {
            use_uniform_crossover = true;
        } else if (arg == "--baseline") {
            output_baseline_edges = true;
        }
    }

    // check mandatory arguments


    if (vertex_filename == "" || edge_filename == "" || train_pos_filename == "" ||
        train_neg_filename == "" || output_dir == "" || str_vertex_events == "" ||
        edge_attributes_type == "")
    {
        std::cerr << "Options --vertices --edges --train_pos --train_neg -o -e -a are mandatory." << std::endl;
        return 1;
    }

    run_process_new(vertex_filename, edge_filename,
                    train_pos_filename, train_neg_filename,
                    output_dir, vertex_events, edge_attributes_type,
                    n_patterns, use_sample, n_sample,
                    test_pos_filename, test_neg_filename,
                    undirected, max_pattern_edges, random_walks, prob_restart,
                    time_unit_primary, time_unit_secondary,
                    evolution_epochs, evolution_subepochs, evaluation_random_walks,
                    use_uniform_crossover, output_baseline_edges, use_seed, seed,
                    use_vertex_attributes, use_simple_init, limit_negative_population,
                    verbose);

    return 0;
}






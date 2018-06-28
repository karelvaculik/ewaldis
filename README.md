EWALDIS
===================

An algorithm for mining discriminative patterns in dynamic graphs.
Accompanying paper was submitted to Discovery Science (DS) 2018 conference.

## Building EWALDIS (in the project directory) in Linux

- Requires **cmake** of version >= 3.9 and c++17 compiler

```sh
cmake .
make
```
## Running the experiments for DS 2018 paper

- Results will be written into results/* directories

```sh
./ewaldis -v --vertices data/dblp1/dblp_vertices.csv --edges data/dblp1/dblp_edges.csv --train_pos data/dblp1/dblp_train_positive_events.csv --train_neg data/dblp1/dblp_train_negative_events.csv -o results/dblp1 -e edges -a nominal -p 20 -n 30 --test_pos data/dblp1/dblp_test_positive_events.csv --test_neg data/dblp1/dblp_test_negative_events.csv -u -m 5 -w 1000 -r 0.3 -pu 1.0 -su 0.5 -ep 25 -se 25 -ew 10 --unicross --baseline -s 1

./ewaldis -v --vertices data/dblp2/dblp_vertices.csv --edges data/dblp2/dblp_edges.csv --train_pos data/dblp2/dblp_train_positive_events.csv --train_neg data/dblp2/dblp_train_negative_events.csv -o results/dblp2 -e edges -a nominal -p 20 -n 30 --test_pos data/dblp2/dblp_test_positive_events.csv --test_neg data/dblp2/dblp_test_negative_events.csv -u -m 5 -w 1000 -r 0.3 -pu 1.0 -su 0.5 -ep 25 -se 25 -ew 10 --unicross --baseline -s 1

./ewaldis -v --vertices data/dblp3/dblp_vertices.csv --edges data/dblp3/dblp_edges.csv --train_pos data/dblp3/dblp_train_positive_events.csv --train_neg data/dblp3/dblp_train_negative_events.csv -o results/dblp3 -e edges -a nominal -p 20 -n 30 --test_pos data/dblp3/dblp_test_positive_events.csv --test_neg data/dblp3/dblp_test_negative_events.csv -u -m 5 -w 1000 -r 0.3 -pu 1.0 -su 0.5 -ep 25 -se 25 -ew 10 --unicross --baseline -s 1

./ewaldis -v --vertices data/enron/enron_vertices.csv --edges data/enron/enron_edges.csv --train_pos data/enron/enron_train_positive_events.csv --train_neg data/enron/enron_train_negative_events.csv -o results/enron -e edges -a nominal -p 10 -n 10 --test_pos data/enron/enron_test_positive_events.csv --test_neg data/enron/enron_test_negative_events.csv -m 5 -w 3000 -r 0.3 -pu 1728000 -su 345600 -ep 25 -se 25 -ew 10 --unicross --baseline -s 1
```

## Description of EWALDIS output files

- For each experiment, there is a directory with the following three files for each discriminative pattern:
- `results_K.txt` with information about the pattern itself:
    - DATA - name of the experiment
    - TOTAL_TIME_PREPARATION - running time of data preparation in seconds (this phase is common for all patterns from one experiment and thus should be counted only once in total time)
    - TOTAL_TIME_MINING - running time of pattern mining in seconds
    - TOTAL_TIME_EVALUATION - running time of pattern evaluation in seconds
    - TIME_UNIT_PRIMARY - primary time unit parameter
    - TIME_UNIT_SECONDARY - secondary time unit parameter
    - RANDOM_WALKS - number of random walks
    - PROB_RESTART - probability of restart
    - MAX_PATTERN_EDGES - maximum number of pattern edges to be found
    - EVOLUATION_EPOCHS - number of epochs in the genetic algorith
    - EVOLUATION_SUBEPOCHS - number of subepochs in the genetic algorith
    - EVALUATION_RANDOM_WALKS - number of random walks used for evaluation
    - PATTERN_EDGES - list of pattern edges in the form of (source_vertex, destination_vertex)
    - PATTERN_SCORES - scores of pattern edges (their fitness from the genetic algorithm)
    - PATTERN_TIMESTAMPS - distribution of relative timestamps for each pattern edge
    - PATTERN_ATTRIBUTES - distribution of attributes for each pattern edge
    - PATTERN_DIRECTIONS - True = edge directed from source to destination, False = the opposite direction; makes sense only for directed graph
    - PATTERN_UNDIRECTED - whether the patten is an undirected graph
    - SAMPLE_EVALUATION_POSITIVE - matching scores on positive events from sample set
    - SAMPLE_EVALUATION_NEGATIVE - matching scores on negative events from sample set
    - TRAIN_EVALUATION_POSITIVE - matching scores on positive events from train set
    - TRAIN_EVALUATION_NEGATIVE - matching scores on negative events from train set
    - TEST_EVALUATION_POSITIVE - matching scores on positive events from test set
    - TEST_EVALUATION_NEGATIVE - matching scores on negative events from test set
- `results_K_positive_population.txt`
    - For each pattern edge, there are summary statistics of fitness score for each epoch
    - For example, if there are 20 epochs, then first 20 rows belong to the first pattern edge, etc.
    - Each row contains 0%, 10%, 50%, 90%, and 100% percentile of the fitness
    - There are exactly MAX_PATTERN_EDGES * EVOLUATION_EPOCHS rows in the file; if the algorithm finished earlier for a given pattern, then there are missing value
- `results_K_negative_population.txt`
    - Similar information as in the previous file, but now for each pattern edge and each epoch, there are summary statistics of negative population fitness for each subepoch
    - Thus, there are MAX_PATTERN_EDGES * EVOLUATION_EPOCHS * EVOLUATION_SUBEPOCHS rows
    - Again, if the computation finishes earlier, there may be missing values
- Moreover, there are two csv files used for creating baseline dataset by R scripts: `baseline_edges_data_train.csv` and `baseline_edges_data_test.csv`

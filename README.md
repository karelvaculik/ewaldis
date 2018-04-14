EWALDIS
===================

An algorithm for mining discriminative patterns in dynamic graphs.
Accompanying paper was submitted to ECML PKDD 2018 conference.

## Building EWALDIS (in the project directory) in Linux

- Requires **cmake** of version >= 3.9 and c++17 compiler

```sh
cmake .
make
```

## Running the experiments for ECML PKDD 2018 paper

- Data files `dblp_edges.csv`, `dblp_vertices.csv`, `enron_edges.csv`, and `enron_vertices.csv` are assumed to exist in `data` directory
- EWALDIS used for pattern extraction (this will create directory `results` with all output files):

```sh
# by selecting specific experiment
./ewaldis dblp1

# or you can run all experiments at once:
./ewaldis dblp1 dblp2 dblp3 dblp4 enron

# or just specific experiments:
./ewaldis enron dblp3 dblp4
```

- Please note that the results may be a little bit different from the ones published in ECML PKDD paper due to randomness in the computation
- To create ARFF files from the EWALDIS results, use `scripts/dataset_preparation*.R` scripts
- In order to run `dataset_preparation_dblp.R` (or `dataset_preparation_enron.R`) correctly,
you have to set working directory and/or paths to results correctly in R environment
- If run correctly, the script will output the following files (similarly for enron):
   - `dblp_ewaldis_train.arff`
   - `dblp_ewaldis_test.arff`
   - `dblp_baseline_train.arff`
   - `dblp_baseline_test.arff`
- These `arff` files can be loaded into [WEKA](https://www.cs.waikato.ac.nz/ml/weka/) for example

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


## Description of ARFF files
- The R scripts create 4 ARFF files
    - `*_ewaldis_train.arff` - for each pattern, there is one attribute with matching score of that pattern on the given event from training set (i.e. there is one row for each event). There are two classes: `pos` and `neg`
    - `*_ewaldis_test.arff` - the same, but now for test events
    - `*_baseline_train.arff` - the instances are the same as for `*_ewaldis_train.arff`, but attributes were created by the baseline method. For each event, the method looked at all past edges adjacent to this event and encoded the presence of these edges by 0/1 values. Specifically, there is one attribute for each such an edge in the form `fA_T`, where `A` represents edge label (such as "nips") and `T` represents edge time relative to the event (e.g. `10` means 10 time units before the event). 
    - `*_baseline_test.arff` - the same, but now for test events

## Further remarks

- Interface for arbitrary input files is going to be added in a near future

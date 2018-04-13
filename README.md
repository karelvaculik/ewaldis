EWALDIS
===================

An algorithm for mining discriminative patterns in dynamic graphs.
Accompanying paper was submitted to ECML PKDD 2018 conference.

## Installation (from the project directory)

- requires cmake of version >= 3.9 and c++17 compiler

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

## Further remarks

- Interface for arbitrary input files is going to be added in a near future

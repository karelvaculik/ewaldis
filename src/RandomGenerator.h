//
// Created by karbal on 6.1.18.
//

#ifndef EWALDIS_CPP_RANDOMGENERATOR_H
#define EWALDIS_CPP_RANDOMGENERATOR_H

#include <algorithm>
#include <vector>
#include <set>
#include <random>
#include "commonutils.h"

class RandomGenerator
{
private:
    std::mt19937 rng;
public:
    RandomGenerator();
    RandomGenerator(unsigned seed);
    int generate_random_int(int min, int max);
    double generate_random_double(double min, double max);
    int generate_random_int_from_distribution(std::vector<double> probs);
    std::pair<int, int> generate_random_int_pair(int max);
    std::vector<int> generate_random_int_vector(int n, int max);

    std::vector<int> stochastic_universal_sampling(int n, std::vector<double> population);

};


#endif //EWALDIS_CPP_RANDOMGENERATOR_H

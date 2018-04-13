//
// Created by karbal on 6.1.18.
//

#include "RandomGenerator.h"

#include <chrono>

//#include <iostream>


RandomGenerator::RandomGenerator() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    rng = std::mt19937(seed);
}

RandomGenerator::RandomGenerator(unsigned seed) {
    rng = std::mt19937(seed);
}


int RandomGenerator::generate_random_int(int min, int max) {
    // randomly generates an integer from interval [min, max-1]
    std::uniform_int_distribution<int> uni(min, max-1);
    return uni(this->rng);
}

double RandomGenerator::generate_random_double(double min, double max)
{
    std::uniform_real_distribution<double> dist(min, max);
    return dist(this->rng);
}

int RandomGenerator::generate_random_int_from_distribution(std::vector<double> probs) {
    // randomly generates an integer from interval [0, probs.size()-1] by using the probs as probabilities
    std::discrete_distribution<> d(probs.begin(), probs.end());
    return d(this->rng);
}


std::pair<int, int> RandomGenerator::generate_random_int_pair(int max)
{
    // generates two nonequal integers from interval [0, max-1]
    std::vector<int> v = generate_integer_sequence(max);
    std::shuffle(v.begin(), v.end(), rng);
    return std::make_pair(v.at(0), v.at(1));
};


std::vector<int> RandomGenerator::generate_random_int_vector(int n, int max)
{
    /*
     * generates n pair-wise nonequal integers from interval [0, max-1]
     * if n > max, then it creates floor(n / max) concatenated shuffles of values [0, max-1] and also appends
     * n - (floor(n / max) * max) pair-wise nonequal integers from interval [0, max-1]
     * examples:
     * generate_random_int_vector(0, 5) -> []
     * generate_random_int_vector(3, 5) -> [4, 2, 3]
     * generate_random_int_vector(7, 5) -> [3, 2, 0, 1, 4, 1, 0]
     * generate_random_int_vector(13, 5) -> [4, 2, 1, 0, 3, 1, 3, 2, 4, 0, 0, 2, 1]
     * generate_random_int_vector(15, 5) -> [0, 4, 1, 2, 3, 2, 0, 3, 1, 4, 3, 2, 0, 4, 1]
    */
    std::vector<int> all_values;
    while (n > max)
    {
        std::vector<int> v = generate_integer_sequence(max);
        std::shuffle(v.begin(), v.end(), rng);
        all_values.insert( all_values.end(), v.begin(), v.end() );
        n -= max;
    }
    if (n > 0)
    {
        std::vector<int> v = generate_integer_sequence(max);
        std::shuffle(v.begin(), v.end(), rng);
        // truncate vector v to maximum size n
        v.resize(std::min(max, n));
        all_values.insert( all_values.end(), v.begin(), v.end() );
    }
    return all_values;
};



std::vector<int> RandomGenerator::stochastic_universal_sampling(int n, std::vector<double> population)
{
    // https://en.wikipedia.org/wiki/Stochastic_universal_sampling
    // inspired by: http://puzzloq.blogspot.cz/2013/03/stochastic-universal-sampling.html
    double f = 0.0;
    for (int i = 0; i < population.size(); ++i)
    {
        f += population.at(i);
    }
    double p = f / n;
    double rand_double = generate_random_double(0.0, 1.0);
    double start = rand_double * p;

    std::set<int> individiuals_set;

    int index = 0;
    double sum = population.at(index);

    for (int j = 0; j < n; ++j) {
        double pointer = start + j * p;

        if (sum >= pointer) {
            individiuals_set.insert(index);
        } else {
            for (++index; index < population.size(); index++) {
                sum += population[index];
                if (sum >= pointer) {
                    individiuals_set.insert(index);
                    break;
                }
            }
        }
    }
    std::vector<int> individuals( individiuals_set.begin(), individiuals_set.end() );
    return individuals;
}




//
// Created by karbal on 29.12.17.
//

#include <cmath>
#include "RandomWalker.h"
#include "commonutils.h"

std::vector<int> generate_integer_sequence(int n)
{
    // generates a vector of numbers 0 to (n-1)
    std::vector<int> elements(n);
    std::generate(elements.begin(), elements.end(), [i = 0] () mutable { return i++; });
    return elements;
}


std::vector<double> extract_quantiles(std::vector<double> & values, std::vector<double> probs)
{
    // for empty input vector return empty vector

    if (values.size() == 0) return std::vector<double>();

    int n = values.size() - 1;

    // sort values
    std::sort(values.begin(), values.end());

    std::vector<double> results;
    for (int i = 0; i < probs.size(); ++i)
    {
        int index1 = floor(n * probs[i]);
        int index2 = ceil(n * probs[i]);
        results.push_back((values[index1] + values[index2]) / 2.0);
    }
    return results;
}
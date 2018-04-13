//
// Created by karel on 25.12.17.
//

#ifndef EWALDIS_CPP_COMMONUTILS_H
#define EWALDIS_CPP_COMMONUTILS_H

#include <algorithm>
#include <iostream>
#include <cstdarg>
#include <vector>

enum AttributeType { NUMERIC, NOMINAL };

typedef double timestamp_t;

namespace constants
{
    const int INVALID_ID = -1;
}

std::vector<int> generate_integer_sequence(int n);
std::vector<double> extract_quantiles(std::vector<double> & values, std::vector<double> probs);

#endif //EWALDIS_CPP_COMMONUTILS_H

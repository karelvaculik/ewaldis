//
// Created by karel on 25.12.17.
//

#include <iostream>
#include "NominalEncoder.h"

double NominalEncoder::get_encoding(std::string value)
{
    if (value_to_encoding.count(value) > 0)
    {
        return (double) value_to_encoding.at(value);
    }
    else
    {
        int new_encoding = value_to_encoding.size();
        value_to_encoding[value] = new_encoding;
        encoding_to_value[new_encoding] = value;
        return (double) new_encoding;
    }
}

std::string NominalEncoder::get_value(double encoding)
{
    int key = (int) encoding;
    return encoding_to_value.at(key);
}
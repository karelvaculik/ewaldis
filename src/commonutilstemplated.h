//
// Created by karel on 22.1.18.
//

#ifndef EWALDIS_CPP_COMMONUTILSPRINT_H
#define EWALDIS_CPP_COMMONUTILSPRINT_H

#include "Edge.h"
#include "Vertex.h"
#include <iostream>
#include <vector>
#include <map>

using namespace std;

template <typename K, typename V> inline bool map_contains(const std::map<K, V> & t, K el);


template <typename K, typename V> bool map_contains(const std::map<K, V> & t, K el);
inline std::string str(Edge * t);
inline std::string str(Edge t);
inline std::string str(Vertex * t);
inline std::string str(Vertex t);
template<typename T>
enable_if_t<is_arithmetic<T>::value, string> str(const T& t);
template<typename T>
enable_if_t<!is_arithmetic<T>::value, string> str(const T& t);
template <typename T1, typename T2> std::string str(const std::pair<T1, T2> & t);
template <typename T> std::string str(const std::vector<T> & t);
template <typename T1, typename T2> std::string str(const std::map<T1, T2> & t);
template <typename First, typename... Rest> std::string str(const First& first, const Rest&... rest);
inline void print();
inline void println();
template <typename T> void print(const T& t);
template <typename T> void println(const T& t);
template <typename First, typename... Rest> void print(const First& first, const Rest&... rest);
template <typename First, typename... Rest> void println(const First& first, const Rest&... rest);


#include <commonutilstemplated.cpp>


#endif //EWALDIS_CPP_COMMONUTILSPRINT_H

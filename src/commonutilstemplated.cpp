//
// Created by karel on 22.1.18.
//


#ifdef EWALDIS_CPP_COMMONUTILSPRINT_H

#include <type_traits>
#include <sstream>
#include "commonutilstemplated.h"


template <typename K, typename V> bool map_contains(const std::map<K, V> & t, K el)
{
    // checks whether the map contains given key:
    return (t.find( el ) != t.end());
}

//inline void print() {
//    cout << endl;
//}
//inline void println() {
//    cout << endl;
//}


//std::string str() {
//    return "";
//}



//inline void print(Edge * t) {
//    print("<e:", t->get_edge_id(), "|", t->get_original_edge_id(), ">");
//}
//
//inline void println(Edge * t) {
//    print(t);
//    cout << endl;
//}

//inline void print(Edge * t) {
//    print("<e:", t->get_edge_id(), "|", t->get_original_edge_id(), "|", t->get_from_vertex_id(), "->",
//        t->get_to_vertex_id(), "|", t->get_attributes());
//    print("|", t->get_timestamp(), ">");
//}
//
//inline void println(Edge * t) {
//    print(t);
//    cout << endl;
//}

//template <typename T> std::string str(const NumericType& t) {
//    return std::to_string(t);
//}
//
//template <typename T> std::string str(const T& t) {
//    return "" + t;
//}



inline std::string str(Edge * t)
{
    return "<e:" + str(t->get_edge_id()) + "|" + str(t->get_original_edge_id()) + "|" + str(t->get_from_vertex_id()) +
           "->" + str(t->get_to_vertex_id()) + "|" + str(t->get_attributes()) + "|" +  str(t->get_timestamp()) + ">" ;
}
inline std::string str(Edge t)
{
    return "<e:" + str(t.get_edge_id()) + "|" + str(t.get_original_edge_id()) + "|" + str(t.get_from_vertex_id()) +
           "." + str(t.get_to_vertex_id()) + "|" + str(t.get_attributes()) + "|" +  str(t.get_timestamp()) + ">" ;
}
inline std::string str(Vertex * t)
{
    return "<v:" + str(t->get_vertex_id()) + "|" + str(t->get_attributes()) + ">" ;
}
inline std::string str(Vertex t)
{
    return "<v:" + str(t.get_vertex_id()) + "|" + str(t.get_attributes()) + ">" ;
}

template<typename T>
enable_if_t<is_arithmetic<T>::value, string> str(const T& t) {
    return std::to_string(t);
//    return std::to_string(t);
}

template<typename T>
enable_if_t<!is_arithmetic<T>::value, string> str(const T& t) {
//    return "" + t;
    return static_cast<ostringstream&>(ostringstream() << t).str();
}


//template <typename T> void print(const T& t) {
//    cout << t;
//}
//template <typename T> void println(const T& t) {
//    cout << t << endl;
//}



//template <typename First, typename... Rest> std::string str(const First& first, const Rest&... rest) {
//    return str(first) + str(rest...);
//}
//
//template <typename First, typename... Rest> void print(const First& first, const Rest&... rest) {
//    cout << first;
//    print(rest...); // recursive call using pack expansion syntax
//}
//
//template <typename First, typename... Rest> void println(const First& first, const Rest&... rest) {
//    cout << first;
//    println(rest...); // recursive call using pack expansion syntax
//}




template <typename T1, typename T2> std::string str(const std::pair<T1, T2> & t) {
    return "(" + str(t.first) + ", " + str(t.second) + ")";
}


//template <typename T1, typename T2> void print_pair(const std::pair<T1, T2> & t)
//{
//    print("(");
//    print(t.first);
//    print(", ");
//    print(t.second);
//    print(")");
//}


//template <typename T1, typename T2> void print(const std::pair<T1, T2> & t) {
//    print_pair(t);
//}
//
//template <typename T1, typename T2> void println(const std::pair<T1, T2> & t) {
//    print_pair(t);
//    println("");
//}


template <typename T> std::string str(const std::vector<T> & t) {
    std::string output = "[";
    for (int i = 0; i < t.size(); ++i) {
        output += str(t[i]);
//        print(t[i]);
        if (i < t.size() - 1)
        {
            output += ", ";
//            cout << ", ";
        }
    }
//    print("]");
    output += "]";
    return output;
}


//template <typename T> void print_vector(const std::vector<T> & t)
//{
//    print("[");
//    for (int i = 0; i < t.size(); ++i) {
//        print(t[i]);
//        if (i < t.size() - 1)
//        {
//            cout << ", ";
//        }
//
//    }
//    print("]");
//}
//
//
//template <typename T> void print(const std::vector<T> & t) {
//    print_vector(t);
//}
//
//template <typename T> void println(const std::vector<T> & t) {
//    print_vector(t);
//    println("");
//}

template <typename T1, typename T2> std::string str(const std::map<T1, T2> & t)
{
//    print("{");
    std::string output = "{";
    int i = 0;
    for (auto & kv : t)
    {
//        print(kv.first);
//        print(": ");
//        print(kv.second);
        output += str(kv.first) + ": " + str(kv.second);
        if (i < t.size() - 1)
        {
//            cout << ", ";
            output += ", ";
        }

        i++;
    }
//    print("}");
    output += "}";
    return output;
}


//template <typename T1, typename T2> void print_map(const std::map<T1, T2> & t)
//{
//    print("{");
//    int i = 0;
//    for (auto & kv : t)
//    {
//        print(kv.first);
//        print(": ");
//        print(kv.second);
//        if (i < t.size() - 1)
//        {
//            cout << ", ";
//        }
//
//        i++;
//    }
//    print("}");
//}
//
//
//
//
//template <typename T1, typename T2> void print(const std::map<T1, T2> & t) {
//    print_map(t);
//}
//
//template <typename T1, typename T2> void println(const std::map<T1, T2> & t) {
//    print_map(t);
//    println("");
//}

//inline std::string str(Edge * t)
//{
//
//    return "<e:" + std::to_string(t->get_edge_id()) + "|" + std::to_string(t->get_original_edge_id()) + ">";
//}



template <typename First, typename... Rest> std::string str(const First& first, const Rest&... rest) {
    return str(first) + str(rest...);
}


inline void print() {
    cout << endl;
}
inline void println() {
    cout << endl;
}

template <typename T> void print(const T& t) {
    cout << str(t);
}
template <typename T> void println(const T& t) {
    cout << str(t) << endl;
}


template <typename First, typename... Rest> void print(const First& first, const Rest&... rest) {
//    cout << first;
//    print(rest...); // recursive call using pack expansion syntax

//    cout << str(first) << str(rest...);
    cout << str(first, rest...);

}

template <typename First, typename... Rest> void println(const First& first, const Rest&... rest) {
//    cout << first;
//    println(rest...); // recursive call using pack expansion syntax

    cout << str(first, rest...) << endl;

}



#endif

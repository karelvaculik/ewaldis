project(src)

SET(HEADER_FILES
        commonutils.h
        Vertex.h
        Edge.h
        DynamicGraph.h
        DynamicGraphExamples.h
        NominalEncoder.h
        RandomWalker.h
        PatternMiner.h
        Pattern.h
        Suitabilities.h
        GeneticPatternExtractor.h
        PatternEdge.h
        OccupiedGraph.h
        RandomGenerator.h
        PatternChecker.h
        commonutilstemplated.h
        )

set(SOURCE_FILES
        commonutils.cpp
        Vertex.cpp
        Edge.cpp
        DynamicGraph.cpp
        DynamicGraphExamples.cpp
        NominalEncoder.cpp
        RandomWalker.cpp
        PatternMiner.cpp
        Pattern.cpp
        Suitabilities.cpp
        GeneticPatternExtractor.cpp
        PatternEdge.cpp
        OccupiedGraph.cpp
        RandomGenerator.cpp
        PatternChecker.cpp
        commonutilstemplated.cpp)

add_library(src STATIC ${SOURCE_FILES} ${HEADER_FILES})
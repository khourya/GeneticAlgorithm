#pragma once

#include <random>

class Generation
{
    int m_size = -1;
};

enum class SelectionType
{
    Choice = 0,  // per Selected Fit Individual:
                 // go pick one other Individual randomly from the entire population
                 // crossover to create 1 child which makes it in the next population

    MaxFit       // Pair the two successively most fit Individuals in fittestIndividuals
                 // Cross over to create two children which make it in the next population
};

enum class MutationType
{
    None = 0,
    Jump,
    Creep
};

struct Options
{
    int seed = 0;                                   // psuedo-random seed
    int nIndividuals = 0;                           // number of individuals in a generation
    int nGenerations = 0;                           // number of generations to evaluate
    double fitnessTolerance = 0.0;                  // convergence tolerance, not currently implemented
    double mutationProbability = 0.0;               // probability of a mutation occuring
    MutationType mutationType = MutationType::Jump; // how mutation operates

    std::string resultsLocation = "";
};

namespace BaseFunctions
{
    bool WeightedCoinFlip(double p);
};
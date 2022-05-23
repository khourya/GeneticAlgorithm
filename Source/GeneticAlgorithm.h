#pragma once

#include <iostream>
#include <cfloat>
#include <vector>

#include "DesignVariable.h"
#include "BaseTypesFunctions.h"
#include "Individuals.h"
// #include "MyMap.h"
// #include "ObjectiveFunction.h"

class GeneticAlgorithm
{
public:
    GeneticAlgorithm();
    GeneticAlgorithm(int populationSize, int maxGenerations, double fitnessTolerance, double mutationProbability, MutationType mutationType);
    void Initialize(std::vector<DesignVariable*> designVariables);

    std::vector<DesignVariable*> GetDesignVariables() { return m_designVariables; }

    void CreateInitialPopulation();
    void EvaluateFitness();
    void SelectionAndCrossover(SelectionType method);

    std::vector<Individual*> GetPopulation() { return m_population; }
    int GetCurrentGeneration() { return m_currentGeneration; }
    int GetCurrentGenerationSize() { return m_populationSize; }
    int GetMaxGenerations() { return m_nMaxGenerations; }

    double GetTotalFitnessCurrentGeneration() { return m_currentPopulationTotalFitness; };
    double GetAverageFitnessCurrentGeneration() { return m_currentPopulationAverageFitness; };
    double GetMaxFitnessCurrentGeneration() { return m_currentPopulationMaxFitness; };
    double GetMinFitnessCurrentGeneration() { return m_currentPopulationMinFitness; };
    // std::vector<Individual*> SelectFittest(Map* fitnessMap);

    Individual* GetRandomIndividual();

    // Individual Creation Routines
    Individual* CreateIndividual();

private:
    // Genetic Algorithm Parameters
    int m_initialPopulationSize = -1;
    int m_populationSize = -1;

    int m_nMaxGenerations = -1;
    int m_currentGeneration = -1;

    double m_fitnessTolerance = DBL_MIN;
    MutationType m_mutationType = MutationType::None;
    double m_mutationProbability = DBL_MIN;

    // Design Variable Parameters
    std::vector<DesignVariable*> m_designVariables;
    int m_nVars = -INT_MIN;
    std::vector<int> m_geneLengths;

    // Fitness tracking
    double m_currentPopulationTotalFitness = DBL_MIN;
    double m_currentPopulationAverageFitness = DBL_MIN;
    double m_currentPopulationMaxFitness = DBL_MIN;
    double m_currentPopulationMinFitness = DBL_MIN;

    std::vector<Individual*> m_population;
};
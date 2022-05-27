#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm()
{
	// we're going to do nothing here
}

void GeneticAlgorithm::SetDesignVariables(std::vector<DesignVariable*> designVariables)
{
	m_nVars = static_cast<int>(designVariables.size());
	m_designVariables = designVariables;
}

void GeneticAlgorithm::SetPerformanceParameters(int populationSize, int maxGenerations, double fitnessTolerance, double mutationProbability, MutationType mutationType)
{
	m_initialPopulationSize = populationSize;
	m_nMaxGenerations = maxGenerations;
	m_fitnessTolerance = m_fitnessTolerance;
	m_mutationProbability = mutationProbability;
	m_mutationType = mutationType;
}


void GeneticAlgorithm::CreateInitialPopulation()
{
	for (int i = 0; i < m_initialPopulationSize; i++)
	{
		// TODO: use shared(?) pointers
		Individual* individual = CreateIndividual();
		m_population.push_back(individual);
	}

	m_populationSize = static_cast<int>(m_population.size());
	m_currentGeneration = 1;
}

Individual* GeneticAlgorithm::CreateIndividual()
{
	Individual* individual = new Individual(m_designVariables);
	return individual;
}
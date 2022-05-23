#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm()
{
	// we're going to do nothing here
}

void GeneticAlgorithm::Initialize(std::vector<DesignVariable*> designVariables)
{
	m_nVars = designVariables.size();
	m_designVariables = designVariables;
}

Individual* GeneticAlgorithm::CreateIndividual()
{
	Individual* individual = new Individual(m_designVariables);
	return individual;
}
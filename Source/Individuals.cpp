#include "Individuals.h"

Individual::Individual()
{
    for (int i = 0; i < m_chromosomalLength; i++)
    {
        int x = rand() % 100;
        if (x > 50)
            m_chromosome[i] = 1;
    }
}

Individual::Individual(double wi, double hi, double Ki)
{
    std::vector<short> gene1 = decimalToBinary(wi, Variable::wi);
    std::vector<short> gene2 = decimalToBinary(hi, Variable::hi);
    std::vector<short> gene3 = decimalToBinary(Ki, Variable::Ki);

    std::copy(gene1.begin(), gene1.end(), m_chromosome.begin());
    std::copy(gene2.begin(), gene2.end(), m_chromosome.begin() + m_geneLength);
    std::copy(gene3.begin(), gene3.end(), m_chromosome.begin() + m_geneLength + m_geneLength);
}

Individual::Individual(std::vector<short> parent1, std::vector<short> parent2, double pMutation = 0, MutationType mutationType = MutationType::None)
{
    int ind = rand() % (m_chromosomalLength + 1);
    std::copy(parent1.begin(), parent1.begin() + ind, m_chromosome.begin());
    std::copy(parent2.begin() + ind, parent2.end(), m_chromosome.begin() + ind);

    // TODO: Strategy pattern for mutation!
    if (mutationType == MutationType::None)
        return;
    else if (mutationType == MutationType::Creep)
    {
        for (int i = 0; i < m_chromosomalLength; i++)
        {
            bool b_mutate = BaseFunctions::WeightedCoinFlip(pMutation);
            if (b_mutate)
                m_chromosome[i] = swap(m_chromosome[i]);
        }
    }
    else if (mutationType == MutationType::Jump)
    {
        bool b_mutate = BaseFunctions::WeightedCoinFlip(pMutation);
        if (b_mutate)
            m_chromosome = swap(m_chromosome);
    }
}

void Individual::PrintChromosome()
{
    for (short elem : m_chromosome)
        std::cout << elem;

    std::cout << std::endl;
}

std::string Individual::GetChromosomeString()
{
    std::stringstream ss;
    for (short bit : m_chromosome)
    {
        ss << bit;
    }

    std::string str = ss.str();
    return str;
}

std::vector<short> Individual::GetChromosomeVector()
{
    return m_chromosome;
}

DesignVariables Individual::GetDesignVariables()
{
    DesignVariables designVariables;

    std::vector<short>::const_iterator begin = m_chromosome.begin();
    std::vector<short>::const_iterator split1 = m_chromosome.begin() + m_geneLength;
    std::vector<short>::const_iterator split2 = m_chromosome.begin() + m_geneLength + m_geneLength;
    std::vector<short>::const_iterator end = m_chromosome.end();

    std::vector<short> gene1(begin, split1);
    std::vector<short> gene2(split1, split2);
    std::vector<short> gene3(split2, end);

    double wi = binaryToDecimal(gene1, Variable::wi);
    double hi = binaryToDecimal(gene2, Variable::hi);
    double Ki = binaryToDecimal(gene3, Variable::Ki);
    double hs = m_hTotal - (2. * hi);
    double ws = (m_wTotal - wi) / 2.;

    designVariables.wi = wi;
    designVariables.ws = ws;
    designVariables.hi = hi;
    designVariables.hs = hs;
    designVariables.Ki = Ki;

    return designVariables;
}

double Individual::binaryToDecimal(std::vector<short> gene, Variable variable)
{
    double varMin = -DBL_MAX;
    double varMax = DBL_MAX;

    if (variable == Variable::wi)
    {
        varMin = m_wiMinValue;
        varMax = m_wiMaxValue;
    }
    else if (variable == Variable::hi)
    {
        varMin = m_hiMinValue;
        varMax = m_hiMaxValue;
    }
    else if (variable == Variable::Ki)
    {
        varMin = m_KiMinValue;
        varMax = m_KiMaxValue;
    }

    double base10 = 0.;
    for (size_t i = 0; i < gene.size(); i++)
    {
        size_t j = gene.size() - 1 - i;
        base10 += gene[j] * pow(2, i);
    }

    base10 = varMin + base10 * (varMax - varMin) / (pow(2., gene.size()) - 1.);

    return base10;
}

std::vector<short> Individual::decimalToBinary(double x, Variable variable)
{
    std::vector<short> gene(m_geneLength, 0);
    double minTemp = -DBL_MAX;
    double maxTemp = DBL_MAX;

    if (variable == Variable::wi)
    {
        minTemp = m_wiMinValue;
        maxTemp = m_wiMaxValue;
    }
    else if (variable == Variable::hi)
    {
        minTemp = m_hiMinValue;
        maxTemp = m_hiMaxValue;
    }
    else if (variable == Variable::Ki)
    {
        minTemp = m_KiMinValue;
        maxTemp = m_KiMaxValue;
    }
    
    for (int i = 0; i < m_geneLength; i++)
    {
        double rngTemp = maxTemp - minTemp;
        double midpoint = minTemp + rngTemp / 2.;

        if (x > midpoint)
        {
            gene[i] = 1;
            maxTemp = maxTemp;
            minTemp = midpoint;
        }
        else
        {
            gene[i] = 0;
            maxTemp = midpoint;
            minTemp = minTemp;
        }
    }

    return gene;
}

short Individual::swap(short bit)
{
    switch (bit)
    {
    case 0:
        return 1;
        break;
    case 1:
        return 0;
        break;
    }

    return 0;
}

std::vector<short> Individual::swap(std::vector<short> chromosome)
{
    std::vector<short> newChromosome(m_chromosomalLength, 0);
    for (int i = 0; i < m_chromosomalLength; i++)
    {
        newChromosome[i] = swap(chromosome[i]);
    }

    return newChromosome;
}
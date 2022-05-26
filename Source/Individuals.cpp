#include "Individuals.h"

/*
Individual::Individual()
{
    for (int i = 0; i < m_chromosomalLength; i++)
    {
        int x = rand() % 100;
        if (x > 50)
            m_chromosome[i] = 1;
    }
}

Individual::Individual(std::vector<double> vars)
{
    for (int i = 0; i < m_nVariables; i++)
    {
        int startCopyIndex = i * m_geneLength;
        std::vector<short> gene = decimalToBinary(vars[i], Variable::var1);
        std::copy(gene.begin(), gene.end(), m_chromosome.begin() + startCopyIndex);
    }
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

*/

Individual::Individual(std::vector<DesignVariable*> designVariables)
{
    m_chromosomalLength = 0;
    m_nVariables = static_cast<int>(designVariables.size());

    for (DesignVariable* designVar : designVariables)
    {
        m_chromosomalLength = m_chromosomalLength + designVar->GetLength();
    }

    m_chromosome.resize(m_chromosomalLength);

    for (int i = 0; i < m_chromosomalLength; i++)
    {
        int x = rand() % 100;
        if (x > 50)
            m_chromosome[i] = 1;
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

DesignVariablesStruct Individual::GetDesignVariables()
{
    DesignVariablesStruct designVariables;

    std::vector<short>::const_iterator begin = m_chromosome.begin();
    std::vector<short>::const_iterator split1 = m_chromosome.begin() + m_geneLength;
    std::vector<short>::const_iterator split2 = m_chromosome.begin() + m_geneLength + m_geneLength;
    std::vector<short>::const_iterator end = m_chromosome.end();

    std::vector<short> gene1(begin, split1);
    std::vector<short> gene2(split1, split2);
    std::vector<short> gene3(split2, end);

    double wi = binaryToDecimal(gene1, Variable::var1);
    double hi = binaryToDecimal(gene2, Variable::var2);
    double Ki = binaryToDecimal(gene3, Variable::var3);
    double hs = m_hTotal - (2. * hi);
    double ws = (m_wTotal - wi) / 2.;

    designVariables.var1 = wi;
    designVariables.var2 = ws;
    designVariables.var3 = hi;
    designVariables.var4 = hs;
    designVariables.var5 = Ki;

    return designVariables;
}

double Individual::binaryToDecimal(std::vector<short> gene, Variable variable)
{
    double varMin = -DBL_MAX;
    double varMax = DBL_MAX;

    if (variable == Variable::var1)
    {
        varMin = m_wiMinValue;
        varMax = m_wiMaxValue;
    }
    else if (variable == Variable::var2)
    {
        varMin = m_hiMinValue;
        varMax = m_hiMaxValue;
    }
    else if (variable == Variable::var3)
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

    if (variable == Variable::var1)
    {
        minTemp = m_wiMinValue;
        maxTemp = m_wiMaxValue;
    }
    else if (variable == Variable::var2)
    {
        minTemp = m_hiMinValue;
        maxTemp = m_hiMaxValue;
    }
    else if (variable == Variable::var3)
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

std::vector<short> Individual::decimalToBinary(double x, double x_min, double x_max)
{
    std::vector<short> gene(m_geneLength, 0);
    double varMin = x_min;
    double varMax = x_max;

    for (int i = 0; i < m_geneLength; i++)
    {
        double rngTemp = x_max - x_min;
        double midpoint = x_min + rngTemp / 2.;

        if (x > midpoint)
        {
            gene[i] = 1;
            x_max = x_max;
            x_min = midpoint;
        }
        else
        {
            gene[i] = 0;
            x_max = midpoint;
            x_min = x_min;
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
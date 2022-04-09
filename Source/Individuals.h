#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include "BaseTypesFunctions.h"

enum class Variable
{
    wi = 0,
    hi,
    Ki
};

struct DesignVariables
{
    double wi = -DBL_MAX;
    double ws = -DBL_MAX;
    double hi = -DBL_MAX;
    double hs = -DBL_MAX;
    double Ki = -DBL_MAX;
};

class Individual
{
public:
    Individual();
    Individual(double wi, double hi, double Ki);
    Individual(std::vector<short> parent1, std::vector<short> parent2, double pMutation, MutationType mutationType);
    

    void PrintChromosome();
    std::string GetChromosomeString();
    std::vector<short> GetChromosomeVector();
    DesignVariables GetDesignVariables();
    // std::vector<short> GetBinaryValue();

    void SetFitness(double fitness) { m_fitness = fitness; }
    double GetFitness() { return m_fitness; }

private:
    int m_nVariables = 3;
    double m_fitness = DBL_MIN;
    const int m_chromosomalLength = 24;
    const int m_geneLength = 8;

    const double m_wiMinValue = 0.16;
    const double m_wiMaxValue = 0.22;

    const double m_hiMinValue = 0.48;
    const double m_hiMaxValue = 0.4975;

    const double m_KiMinValue = 0.007;
    const double m_KiMaxValue = 0.1;

    const double m_wTotal = 0.24;
    const double m_hTotal = 1.;

    std::vector<short> m_chromosome = std::vector<short>(m_chromosomalLength, 0);

    std::vector<short> decimalToBinary(double x, Variable variable);
    double binaryToDecimal(std::vector<short>, Variable variable);
    short swap(short bit);
    std::vector<short> swap(std::vector<short> chromosome);
};
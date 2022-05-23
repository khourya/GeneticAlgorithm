#pragma once

#include <iostream>
#include <sstream>
#include <vector>
#include "BaseTypesFunctions.h"
#include "DesignVariable.h"

enum class Variable
{
    var1 = 0,
    var2,
    var3
};

struct DesignVariablesStruct
{
    double var1 = -DBL_MAX;
    double var2 = -DBL_MAX;
    double var3 = -DBL_MAX;
    double var4 = -DBL_MAX;
    double var5 = -DBL_MAX;
};

class Individual
{
public:
    Individual();
    Individual(std::vector<DesignVariable*> designVariables);
    Individual(std::vector<double> vars);
    Individual(std::vector<short> parent1, std::vector<short> parent2, double pMutation, MutationType mutationType);

    void PrintChromosome();
    std::string GetChromosomeString();
    std::vector<short> GetChromosomeVector();
    DesignVariablesStruct GetDesignVariables();
    // std::vector<short> GetBinaryValue();

    void SetFitness(double fitness) { m_fitness = fitness; }
    double GetFitness() { return m_fitness; }

private:
    int m_nVariables = -INT_MIN;
    double m_fitness = DBL_MIN;
    int m_chromosomalLength = -INT_MIN;
    int m_geneLength = -INT_MIN;

    double m_wiMinValue = 0.16;
    double m_wiMaxValue = 0.22;

    double m_hiMinValue = 0.48;
    double m_hiMaxValue = 0.4975;

    double m_KiMinValue = 0.007;
    double m_KiMaxValue = 0.1;

    double m_wTotal = 0.24;
    double m_hTotal = 1.;

    std::vector<short> m_chromosome; // = std::vector<short>(m_chromosomalLength, 0);

    std::vector<short> decimalToBinary(double x, Variable variable);
    std::vector<short> decimalToBinary(double x, double x_min, double x_max);
    double binaryToDecimal(std::vector<short>, Variable variable);
    short swap(short bit);
    std::vector<short> swap(std::vector<short> chromosome);
};
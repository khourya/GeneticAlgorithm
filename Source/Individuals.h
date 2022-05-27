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

class Individual
{
public:
    // Individual();
    Individual(std::vector<DesignVariable*> designVariables);

    // Individual(std::vector<double> vars);
    Individual(std::vector<DesignVariable*> designVariables, std::vector<double> vars);


    // Individual(std::vector<short> parent1, std::vector<short> parent2, double pMutation, MutationType mutationType);

    void PrintChromosome();
    std::string GetChromosomeString();
    std::vector<short> GetChromosomeVector();
    std::vector<DesignVariable> GetDesignVariables();
    // std::vector<short> GetBinaryValue();

    void SetFitness(double fitness) { m_fitness = fitness; }
    double GetFitness() { return m_fitness; }

private:
    int m_nVariables = -INT_MIN;
    double m_fitness = DBL_MIN;
    int m_chromosomalLength = -INT_MIN;
    int m_geneLength = -INT_MIN;

    std::vector<short> m_chromosome; // = std::vector<short>(m_chromosomalLength, 0);
    std::vector<double> m_variableValues;

    std::vector<short> decimalToBinary(double x, DesignVariable* variable);
    std::vector<short> decimalToBinary(double x, double x_min, double x_max);
    double binaryToDecimal(std::vector<short>, DesignVariable* variable);
    short swap(short bit);
    std::vector<short> swap(std::vector<short> chromosome);
};
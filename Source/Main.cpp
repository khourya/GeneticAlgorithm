// Main.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "GeneticAlgorithm.h"
#include "Individuals.h"

int main()
{
    srand(3);
    std::cout << "Hello World!\n";
   
    //                 = new DesignVariable(length, min, max);
    DesignVariable* D1 = new DesignVariable(2, 0.0, 10.0);
    DesignVariable* D2 = new DesignVariable(2, 0.0, 10.0);
    DesignVariable* D3 = new DesignVariable(2, 0.0, 10.0);
    std::vector<DesignVariable*> designVariables = { D1, D2, D3 };

    GeneticAlgorithm* GA = new GeneticAlgorithm();
    GA->Initialize(designVariables);
    Individual* indy = GA->CreateIndividual();
    indy->PrintChromosome();

    std::cout << "check stuff, stud.";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

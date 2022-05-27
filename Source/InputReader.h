#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "BaseTypesFunctions.h"


Options ReadOptions()
{
	std::string line;
	std::ifstream inputFile;
	Options GAOptions;
	int tempint = 0;

	inputFile.open("GAOptions.inp");

	if (!inputFile)
	{
		std::cout << "Unable to open 'options' file." << std::endl;
	}

	std::getline(inputFile, line);    // Read Results Location Path Header
	std::getline(inputFile, line);    // Read Results Location Path
	GAOptions.resultsLocation = line; // setting Results Path

	std::getline(inputFile, line); // skip blank line
	std::getline(inputFile, line); // Read Seed Header
	inputFile >> GAOptions.seed;   // Get seed

	std::getline(inputFile, line);       // skip blank line
	std::getline(inputFile, line);       // skip blank line
	std::getline(inputFile, line);       // Read nGenerations Header
	inputFile >> GAOptions.nGenerations; // Get nGenerations

	std::getline(inputFile, line);       // skip blank line
	std::getline(inputFile, line);       // skip blank line
	std::getline(inputFile, line);       // Read nIndividuals Header
	inputFile >> GAOptions.nIndividuals; // Get nIndividuals

	std::getline(inputFile, line);           // skip blank line
	std::getline(inputFile, line);           // skip blank line
	std::getline(inputFile, line);           // Read Fitness Tolerance Header
	inputFile >> GAOptions.fitnessTolerance; // Get Fitness Tolerance

	std::getline(inputFile, line);              // skip blank line
	std::getline(inputFile, line);              // skip blank line
	std::getline(inputFile, line);              // Read Mutation Probability Header
	inputFile >> GAOptions.mutationProbability; // Get Mutation Probability

	std::getline(inputFile, line); // skip blank line
	std::getline(inputFile, line); // skip blank line
	std::getline(inputFile, line); // Read Mutation Type Header
	inputFile >> tempint;          // Get Mutation Type
	MutationType mutationType = static_cast<MutationType>(tempint);
	GAOptions.mutationType = mutationType;

	inputFile.close();

	return GAOptions;
}
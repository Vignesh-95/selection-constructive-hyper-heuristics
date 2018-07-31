// Some imports might not be used - try eliminate them
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <limits.h>
#include <string>
#include <numeric>

using namespace std;

// One-Dimensional Bin-Packing Problem
struct InputData
{
	int numberOfItems, binCapacity;
	vector<int> items;
};

class Chromosome
{
	public:
		/* Different Representations of Chromosome:
		* Repetition vs Non-Repetition of Heuristics in Chromosome
		* Fixed vs Variable Vs N-times Size Chromosome
		*/	
		vector<int> heuristics;
		// Each element corresponds to a different Training or Test Instance
		vector<int> numberOfBins;
		vector<double> fitnesses;
		// Initialised by GA class.
		// 3 Training and Test Examples - One Per Level Of Difficulty
		static int heuristicsLength;
		static vector<InputData>* input;
	
		/* Characteristics of the Problem:
		* Must be determined during the Packing - while applying different heuristics
		* Stored in these data structures:
		* Length Corresponds to number of bins for particular training/test instance
		*/
		vector<int> binCapacities;
		vector<int> numItemsPerBin;
		
		Chromosome()
		{
			// Random Initialisation of Fixed Length Non-repeating Heuristic Combinations
			for (int count = 0; count < heuristicsLength; count++)
			{
				heuristics.push_back(count);
			}
			srand (time(NULL));
			std::random_shuffle (heuristics.begin(), heuristics.end());
			
			/* Initialise Data Sturctures for training/test instances
			* Calculate Initial Fitnesses By Constructing Solutions using heuristics
			*/
			for (int count = 0; count < 3; count++)
			{
				fitnesses.push_back(-1);
				numberOfBins.push_back(0);
				applyHeuristics(count);
			}	
			for (int count = 3; count < 6; count++)
			{
				fitnesses.push_back(-1);
				numberOfBins.push_back(0);
			}
		}
		
		// Used When Creating Offspring
		void nonInitialisation(vector<int> h)
		{
			for (int count = 0; count < heuristicsLength; count++)
			{
				heuristics.push_back(h[count]);
			}
			
			for (int count = 0; count < 3; count++)
			{
				fitnesses.push_back(-1);
				numberOfBins.push_back(0);
				applyHeuristics(count);
			}

			for (int count = 3; count < 6; count++)
			{
				fitnesses.push_back(-1);
				numberOfBins.push_back(0);
			}			
		}
		
		// Test/Training Instance Number
		void applyHeuristics(int num)
		{
			// Initially a Single Bin - At least one is needed
			numberOfBins[num] = 1;
			binCapacities.clear();
			binCapacities.push_back((*input)[num].binCapacity);
			numItemsPerBin.clear();
			numItemsPerBin.push_back(0);
			// Start from left to right in heuristic combination
			int nextHeuristic = 0;
			
			// Traverse Circular through the heuristic combination till items are exhausted.
			for (int count = 0; count < (*input)[num].numberOfItems; count++)
			{
				/* Pass the training/test instance number and the number of the current item 
				* being packed
				*/
				switch(nextHeuristic)
				{
					case 0: 
						calculateFirstFitFitness(num, count);
						break;
					case 1: 
						calculateBestFitFitness(num, count);
						break;
					case 2: 
						calculateNextFitFitness(num, count);
						break;
					case 3: 
						calculateWorstFitFitness(num, count);
						break;
				}
				
				// Circular Motion
				nextHeuristic = (nextHeuristic + 1) % heuristicsLength;
			}
			
			/* We are done packing
			* Calculate the fitness of the instance using the characteristics of the problem
			*/
			fitnessFunction(num);
		}
		
		/* First-fit (count) - The item to be placed is allocated to the first
		* bin that has sufficient space for it to fit into.
		*/
		void calculateFirstFitFitness(int num, int count)
		{
			bool didFit = false;
			// Cycle through created bins
			for (int b = 0; b < numberOfBins[num]; b++)
			{
				// Can Fit
				if (binCapacities[b]  >= (*input)[num].items[count])
				{
					binCapacities[b] -= (*input)[num].items[count];
					numItemsPerBin[b] += 1;
					didFit = true;
					break;
				}
			}
			
			if (didFit == false)
			{
				binCapacities.push_back((*input)[num].binCapacity - (*input)[num].items[count]);
				numItemsPerBin.push_back(1);
				numberOfBins[num]++;
			}
		}
		
		/* Best-fit (count) - The item is placed into the bin that will have
		* the smallest residual space once the item is placed into it.
		*/
		void calculateBestFitFitness(int num, int count)
		{
			bool didFit = false;
			int minResidue = INT_MAX;
			int minResidueBin = -1;
			for (int b = 0; b < numberOfBins[num]; b++)
			{
				int residue = binCapacities[b]  - (*input)[num].items[count];
				if (residue >= 0 && residue < minResidue)
				{
					minResidue = residue;
					minResidueBin = b;
					// At least one Bin had space i.e. residue >= 0 
					didFit = true;
				}
			}
			
			if (didFit == true)
			{
				binCapacities[minResidueBin] -= (*input)[num].items[count];
				numItemsPerBin[minResidueBin] += 1;
			}
			else 
			{
				binCapacities.push_back((*input)[num].binCapacity - (*input)[num].items[count]);
				numItemsPerBin.push_back(1);
				numberOfBins[num]++;
			}
		}
		
		/* Next-fit (count) - If the item does not fit into the current bin,
		* i.e. the last bin to be created, it is placed in a new bin.
		*/
		void calculateNextFitFitness(int num, int count)
		{
			bool didFit = false;
			if (binCapacities[numberOfBins[num] - 1]  >= (*input)[num].items[count])
			{
				binCapacities[numberOfBins[num] - 1] -= (*input)[num].items[count];
				numItemsPerBin[numberOfBins[num] - 1] += 1;
				didFit = true;
			}
			
			if (didFit == false) 
			{
				binCapacities.push_back((*input)[num].binCapacity - (*input)[num].items[count]);
				numItemsPerBin.push_back(1);
				numberOfBins[num]++;
			}
		}
		
		/* Worst-fit (count) - The item to be placed is allocated to the bin
		* with the largest residual space once the item is placed into it.
		*/
		void calculateWorstFitFitness(int num, int count)
		{
			bool didFit = false;
			int maxResidue = INT_MIN;
			int maxResidueBin = -1;
			for (int b = 0; b < numberOfBins[num]; b++)
			{
				int residue = binCapacities[b]  - (*input)[num].items[count];
				if (residue >= 0 && residue > maxResidue)
				{
					maxResidue = residue;
					maxResidueBin = b;
					didFit = true;
				}
			}
			
			if (didFit == true)
			{
				binCapacities[maxResidueBin] -= (*input)[num].items[count];
				numItemsPerBin[maxResidueBin] += 1;
			}
			else 
			{
				binCapacities.push_back((*input)[num].binCapacity - (*input)[num].items[count]);
				numItemsPerBin.push_back(1);
				numberOfBins[num]++;
			}
		}
		
		void fitnessFunction(int num)
		{
			double sum = 0.0;
			for (int count = 0; count < numberOfBins[num]; count++)
			{
				double calc = (*input)[num].binCapacity - binCapacities[count];
				calc /= (*input)[num].binCapacity;
				calc *= calc;
				
				sum += calc;
			}
			
			fitnesses[num] = sum / numberOfBins[num];
		}
		
		double fitnessOverInputs(bool train)
		{
			double fit = 0.0;
			int start = 0; 
			int end = 3;
			
			if (train != true)
			{
				start = 3;
				end = 6;
			}
			
			for (int count = start; count < end; count++)
			{
				fit += fitnesses[count];
			}
			
			return fit/3.0;
		}
		
		static void setHeuristicsLength(int hLength)
		{
			heuristicsLength = hLength;
		}
		
		static void setInput(vector<InputData>* in)
		{
			input = in;
		}
};

int Chromosome::heuristicsLength;
vector<InputData>* Chromosome::input;

class GA
{
	private:
		vector<Chromosome> population;
		int populationSize, numberOfGenerations, tournamentSize, heuristicsLength;
		double crossoverProbability, mutationProbability;
		// Objective Value For Each Training and Test Instance
		vector<int>* nBins;
		// Best Chromosome Fitness Value
		Chromosome* best;
	
	public:
		GA(int pSize, int ng, int tSize, double cP, double mP, int hl, vector<InputData>* input)
		{
			populationSize = pSize;
			numberOfGenerations = ng;
			tournamentSize = tSize;
			crossoverProbability = cP;
			mutationProbability = mP;
			heuristicsLength = hl;
			
			for (int counter = 0; counter < populationSize; counter++)
			{
				Chromosome::setHeuristicsLength(hl);
				Chromosome::setInput(input);
				population.push_back(Chromosome());
			}
		}
		
		vector<int>* getNBins()
		{
			return nBins;
		}
		
		Chromosome* getBest()
		{
			return best;
		}
		
		static bool compare(Chromosome a, Chromosome b)
		{
			return (a.fitnessOverInputs(true) < b.fitnessOverInputs(true));
		}
		
		int evolve()
		{
			int generationCounter = 0;
			double fitness;
			std::vector<Chromosome>::iterator result;
			while(generationCounter < numberOfGenerations)
			{
				int parents [2];
				selectParents(parents);
				Chromosome* offspring;
				crossover(parents, offspring);
				mutate(offspring);
				insert(parents, offspring);
				generationCounter++;
			}
			
			/*for (int counter = 0; counter < populationSize; counter++)
			{
				for (int count = 0; count < 3; count++)
				{
					population[counter].applyHeuristics(count);
				}
			}*/
			
			result = std::min_element(population.begin(), population.end(), compare);
			best = &(*result);
			fitness = best->fitnessOverInputs(true);
			nBins = &(best->numberOfBins);
			
			return fitness; 
		}
		
		void selectParents(int p [2])
		{
			Chromosome* tournament[tournamentSize];
			vector<int> values;
			int valueLength = 0;
			
			for (int counter = 0; counter < tournamentSize; counter++)
			{
				srand (time(NULL));
				int random = rand() % populationSize;
				
				while (valueLength != 0 && std::find (values.begin(), values.end(), random) != values.end())
				{
					random = rand() % populationSize;
				}
				
				values.push_back(random);
				valueLength++;
				
				tournament[counter] = &population[random];
			}
			
			int parent;
			int parentCounter;
			double min = (double) INT_MAX;
			
			
			/*for (int counter = 0; counter < populationSize; counter++)
			{
				for (int count = 0; count < 3; count++)
				{
					population[counter].applyHeuristics(count);
				}
			}*/
			
			for (int counter = 0; counter < tournamentSize; counter++)
			{
				double fit = tournament[counter]->fitnessOverInputs(true);
				if (fit < min)
				{
					min = fit;
					parent = values[counter];
					parentCounter = counter;
				}
			}
			
			p[0] = parent;
			
			int parent2;
			min = (double) INT_MAX;
			for (int counter = 0; counter < tournamentSize; counter++)
			{	if (counter != parentCounter)
				{
					double fit = tournament[counter]->fitnessOverInputs(true);
					if (fit < min)
					{
						min = fit;
						parent2 = values[counter];
					}
				}	
			}
			p[1] = parent2;
		}
		
		void crossover(int p [2], Chromosome*& offspring)
		{
			int random = rand() % Chromosome::heuristicsLength;
			vector<int> temp1;
			vector<int> temp2;
			for (int it = 0; it < Chromosome::heuristicsLength; it++)
			{
				if (it <= random)
				{	
					temp1.push_back(population[p[0]].heuristics[it]);
					temp2.push_back(population[p[1]].heuristics[it]);
				}
				else
				{
					temp1.push_back(population[p[1]].heuristics[it]);
					temp2.push_back(population[p[0]].heuristics[it]);
				}
			}	
			
			Chromosome* o = new Chromosome;
			o->nonInitialisation(temp1);
			Chromosome* o2 = new Chromosome;
			o->nonInitialisation(temp2);
			
			for (int count = 0; count < 3; count++)
			{
				o->applyHeuristics(count);
				o2->applyHeuristics(count);
			}
			
			if (o->fitnessOverInputs(true) < o2->fitnessOverInputs(true))
			{
				offspring = o;
			}
			else
			{
				offspring = o2;
			}
			
		}
		
		void mutate(Chromosome* offspring)
		{
			srand (time(NULL));
			//if ((rand() % 10000) >= (mutationProbability * 10000))
			//{
			//	return;
			//}
			int random1 = rand() % Chromosome::heuristicsLength;
			int random2 = rand() % Chromosome::heuristicsLength;
			
			offspring->heuristics[random1] = random2;
			
			for (int count = 0; count < 3; count++)
			{
				offspring->applyHeuristics(count);
			}
		}
		
		void insert(int parents[2], Chromosome* offspring)
		{
			double fitnessP1 = population[parents[0]].fitnessOverInputs(true);
			double fitnessP2 = population[parents[1]].fitnessOverInputs(true);
			double fitnessOffspring = offspring->fitnessOverInputs(true);
			
			if (fitnessP1 > fitnessP2)
			{
				if (fitnessP1 > fitnessOffspring)
				{
					population.erase(population.begin() + parents[0]);
					population.push_back(*offspring);
				}
			}
			else
			{
				if (fitnessP2 > fitnessOffspring)
				{
					population.erase(population.begin() + parents[1]);
					population.push_back(*offspring);
				}
			}
		}
};

int main ()
{
	ifstream inFile;
	string fileName;
	int counter = 0;
	vector<InputData> input;
	int x;

	for (int start = 0; start < 6; start++)
	{
		cin >> fileName;
		inFile.open(fileName.c_str());
		input.push_back(InputData());
		
		if (!inFile) 
		{
			cout << "Unable to open file.";
			exit(1);   // call system to stop
		}
		
		while (inFile >> x) 
		{
			if (counter == 0)
			{
				input[start].numberOfItems = x;
			}
			else if (counter == 1)
			{
				input[start].binCapacity = x;
			}
			else
			{
				input[start].items.push_back(x);
			}
			counter++;
		}
		inFile.close();
		counter = 0;
		x = 0;
	}

	// Adjust Parameters
	// Add More Parameters
	GA ga(50, 3, 3, 1, 0.80, 4, &input); 
	
	double fitnessTraining = ga.evolve();
	cout << "\n\n";
	cout << "Fitness Training: " << fitnessTraining << "\n";
	vector<int>* numBinsTraining = ga.getNBins();
	cout << "Number of Bins Training: ";
	for (int q = 0; q < 3; q++)
	{
		cout << (*numBinsTraining)[q] << " ";
	}
	cout << "\n";
	Chromosome* trainingWinner = ga.getBest();
	
	trainingWinner->applyHeuristics(3);
	trainingWinner->applyHeuristics(4);
	trainingWinner->applyHeuristics(5);
	double fitnessTesting = trainingWinner->fitnessOverInputs(false);
	cout << "Fitness Testing: " << fitnessTesting << "\n";
	vector<int>* numBinsTesting = ga.getNBins();
	cout << "Number of Bins Testing: ";
	for (int q = 3; q < 6; q++)
	{
		cout << (*numBinsTesting)[q] << " ";
	}
	cout << "\n\n\n";
	
	return 0;
}
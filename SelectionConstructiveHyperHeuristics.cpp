// Some imports might not be used
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

struct InputData
{
	int numberOfItems, binCapacity;
	vector<int> items;
};

class Chromosome
{
	public:
		vector<int> heuristics;		
		int numberOfBins;
		vector<double> fitnesses;
	
		// Characteristics of Problem must be determined after the Packing
		// Stored in these data structures:
		vector<int> binCapacities;
		vector<int> numItemsPerBin;
		
		// 3 training examples, one from each level.
		// IN GA
		
		static int heuristicsLength;
		static vector<InputData> input;
	
		Chromosome()
		{
			for (int count = 0; count < heuristicsLength; count++)
			{
				heuristics.push_back(count);
			}
			srand (time(NULL));
			std::random_shuffle (heuristics.begin(), heuristics.end());
			
			for (int count = 0; count < 3; count++)
			{
				fitness.push_back(-1);
				applyHeuristics(count);
			}	
		}
		
		void nonInitialisation(vector<int> h)
		{
			for (int count = 0; count < heuristicsLength; count++)
			{
				heuristics.push_back(h[count]);
			}
			
			for (int count = 0; count < 3; count++)
			{
				fitness.push_back(-1);
				applyHeuristics(count);
			}	
		}
		
		void applyHeuristics(int num)
		{
			int numberOfBins = 1;
			binCapacities.clear();
			binCapacities.push_back((*input)[num].binCapacity);
			numItemsPerBin.clear();
			numItemsPerBin.push_back(0);
			int nextHeuristic = 0;
			
			for (int count = 0; count < (*input)[num].numberOfItems; count++)
			{
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
				
				nextHeuristic = (nextHeuristic + 1) % nextHeuristic;
			}
			
			fitnesses[num] = fitnessFunction(num);
		}
		
		void calculateFirstFitFitness(int num, int count)
		{
			bool didFit = false;
			for (int b = 0; b < numberOfBins; b++)
			{
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
				numberOfBins++;
			}
		}
		
		void calculateBestFitFitness(int num, int count)
		{
			bool didFit = false;
			int minResidue = INT_MAX;
			int minResidueBin = -1;
			for (int b = 0; b < numberOfBins; b++)
			{
				int residue = binCapacities[b]  - (*input)[num].items[count];
				if (residue >= 0 && residue < minResidue)
				{
					minResidue = residue;
					minResidueBin = b;
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
				numberOfBins++;
			}
		}
		
		void calculateNextFitFitness(int num, int count)
		{
			bool didFit = false;
			if (binCapacities[numberOfBins - 1]  >= (*input)[num].items[count])
			{
				binCapacities[numberOfBins - 1] -= (*input)[num].items[count];
				numItemsPerBin[numberOfBins - 1] += 1;
				didFit = true;
			}
			
			if (didFit == false) 
			{
				binCapacities.push_back((*input)[num].binCapacity - (*input)[num].items[count]);
				numItemsPerBin.push_back(1);
				numberOfBins++;
			}
		}
		
		void calculateNextFitFitness(int num, int count)
		{
				bool didFit = false;
			int maxResidue = INT_MIN;
			int maxResidueBin = -1;
			for (int b = 0; b < numberOfBins; b++)
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
				numberOfBins++;
			}
		}
		
		double fitnessFunction(int num)
		{
			double sum = 0.0
			for (int count = 0; count < numberOfBins; count++)
			{
				double calc = binCapacities[count];
				calc /= (*input)[num].binCapacity);
				calc *= calc;
				
				sum += calc;
			}
			
			return sum / numberOfBins;
		}
		
		double fitnessOverInputs()
		{
			double fit = 0.0;
			for (int count = 0; count < 3; count++)
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
		
		void printChromosome()
		{
			cout << "\n-----------------------------------------------------------------------\n";
			cout << "\n\nPermutation: ";
			for (int count = 0; count < permutationLength; count++)
			{
				cout << permutation[count] << " ";
			}
			
			cout << "\n\nBin Capacities: ";
			for (int count = 0; count < numberOfBins; count++)
			{
				cout << binCapacities[count] << " ";
			}
			
			cout << "\n\nNumber of Bins: " << numberOfBins;
			cout << "\n-----------------------------------------------------------------------";
		}
		
};

int Chromosome::heuristicsLength;
vector<InputData> Chromosome::input;

class GA
{
	private:
		vector<Chromosome> population;
		int populationSize, numberOfGenerations, tournamentSize;
		double crossoverProbability, mutationProbability;
		int nBins;
	
	public:
		GA(int pSize, int ng, int tSize, double cP, double mP, int hl, vector<InputData>* input)
		{
			populationSize = pSize;
			numberOfGenerations = ng;
			tournamentSize = tSize;
			crossoverProbability = cP;
			mutationProbability = mP;
			
			for (int counter = 0; counter < populationSize; counter++)
			{
				Chromosome::setHeuristicsLength(hl);
				Chromosome::setInput(input);
				population.push_back(Chromosome());
			}
		}
		
		int getNBins()
		{
			return nBins;
		}
		
		static bool compare(Chromosome a, Chromosome b)
		{
			return (a.fitnessOverInputs() < b.fitnessOverInputs());
		}
		
		int evolve()
		{
			int generationCounter = 0;
			double fitness;
			Chromosome* best;
			std::vector<Chromosome>::iterator result;
			while(generationCounter < numberOfGenerations)
			{
				int parents [2];
				selectParents(parents);
				Chromosome* offspring;
				crossover(parents, offspring);
				mutate(offspring);
				insert(parents, offspring);
				//insertDeletingWorst(offspring);
				//best->printChromosome();
				generationCounter++;
			}
			
			result = std::min_element(population.begin(), population.end(), compare);
			best = &(*result);
			fitness = best->fitnessOverInputs();
			nBins = best->getNumberOfBins();
			
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
			for (int counter = 0; counter < tournamentSize; counter++)
			{
				double fit = tournament[counter]->fitnessOverInputs();
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
					double fit = tournament[counter]->fitnessOverInputs();
					if (fit < min)
					{
						min = fit;
						parent2 = values[counter];
					}
				}	
			}
			p[1] = parent2;
		}
		
		void crossover(int p [2], Chromosome* offspring)
		{
			vector<int> off;
			vector<int> seen;
			vector< vector<int> > indicesSeen;
			
			int bigCounter = 0;
			int p1 = p[0];
			int p2 = p[1];
			for (int counter = 0; counter < chromosomeLength; counter++)
			{
				int temp = population[p1].permutationIndex(counter);
				off.push_back(temp);
				vector<int>::iterator it = std::find(seen.begin(), seen.end(), temp);
				if (it != seen.end())
				{
					indicesSeen[it - seen.begin()].push_back(bigCounter);
				}
				else
				{
					seen.push_back(temp);
					vector<int>* t = new vector<int>;
					t->push_back(bigCounter);
					indicesSeen.push_back(*t);
				}
				bigCounter++;
				temp = population[p2].permutationIndex(counter);
				off.push_back(temp);
				it = std::find(seen.begin(), seen.end(), temp);
				if (it != seen.end())
				{
					indicesSeen[it - seen.begin()].push_back(bigCounter);
				}
				else
				{
					seen.push_back(temp);
					vector<int>* t = new vector<int>;
					t->push_back(bigCounter);
					indicesSeen.push_back(*t);
				}
				bigCounter++;
			}
			
			int countUniqueNumbers = seen.size();
			for (int i = 0; i < countUniqueNumbers; i++)
			{
				int innerVectorSize = indicesSeen[i].size()/2;
				for (int j = 0; j < innerVectorSize;  j++)
				{
					off[indicesSeen[i][j]] = INT_MAX;
				}
			}
			
			vector<int> finalOffSpring;
			for (int i = 0; i < bigCounter; i++)
			{
				if (off[i] != INT_MAX)
				{
					finalOffSpring.push_back(off[i]);
				}
			}
			
			Chromosome* o = new Chromosome;
			o->nonInitialisation(finalOffSpring);
			offspring = o;
		}
		
		void mutate(Chromosome* offspring)
		{
			srand (time(NULL));
			//if ((rand() % 10000) >= (mutationProbability * 10000))
			//{
			//	return;
			//}
			
			int random1 = rand() % chromosomeLength;
			int random2 = rand() % chromosomeLength;
			
			int temp = offspring->permutationIndex(random1);
			offspring->setPermutationIndex(random1, offspring->permutationIndex(random2));
			offspring->setPermutationIndex(random2, temp);
		}
		
		void insert(int parents[2], Chromosome* offspring)
		{
			double fitnessP1 = population[parents[0]].fitnessOverInputs();
			double fitnessP2 = population[parents[1]].fitnessOverInputs();
			double fitnessOffspring = offspring->fitnessOverInputs();
			
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
		
		void insertDeletingWorst(Chromosome* offspring)
		{
			std::vector<Chromosome>::iterator result;
			result = std::max_element(population.begin(), population.end(), compare);
			double fitnessP1 = (*result).fitnessOverInputs();
			double fitnessOffspring = offspring->fitnessOverInputs();		
		
			if (fitnessP1 > fitnessOffspring)
			{
				population.erase(result);
				population.push_back(*offspring);
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
	GA ga(100, 100, 3, 1, 0.10, 4, &input); 
	
	// double fitness = ga.evolve();
	// int numBins = ga.getNBins();
	
	return 0;
}
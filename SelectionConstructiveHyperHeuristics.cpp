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

class Chromosome
{
	private:
		vector<int> permutation;
		
		int numberOfBins;
	
		// Characteristics of Problem must be determined after the Packing
		// Stored in these data structures:
		// USE GETTERS
		vector<int> binCapacities;
		vector<int> numItemsPerBin;
		// IN GA
		static int binCapacity;
		static int permutationLength;
	
	public:
		Chromosome()
		{
			
		}
	
		Chromosome(vector<int> items)
		{
			srand (time(NULL));
			std::random_shuffle (items.begin(), items.end());
			for (int count = 0; count < permutationLength; count++)
			{
				permutation.push_back(items[count]);
			}
			calculateFirstFitFitness();
		}
		
		void nonInitialisation(vector<int> items)
		{
			for (int count = 0; count < permutationLength; count++)
			{
				permutation.push_back(items[count]);
			}
			calculateFirstFitFitness();
		}
		
		void calculateFirstFitFitness()
		{
			numberOfBins = 1;
			binCapacities.clear();
			binCapacities.push_back(binCapacity);
			numItemsPerBin.clear();
			numItemsPerBin.push_back(0);
			for (int count = 0; count < permutationLength; count++)
			{
				bool didFit = false;
				for (int b = 0; b < numberOfBins; b++)
				{
					if (binCapacities[b]  >= permutation[count])
					{
						binCapacities[b] -= permutation[count];
						numItemsPerBin[b] += 1;
						didFit = true;
						break;
					}
				}
				if (didFit == false)
				{
					binCapacities.push_back(binCapacity - permutation[count]);
					numItemsPerBin.push_back(1);
					numberOfBins++;
				}
			}
		}
		
		/* TODO: ADD FITNESS FUNCTION FROM GP.
		* RETURN TYPE.
		* PARAMETERS (TERMINAL SET) - IN THIS CLASS
		* BODY (FUNCTION AND TERMINAL SET)
		*/
		double fitnessFunction()
		{
			return 0.0;
		}
		
		static void setBinCapacity(int bCap)
		{
			binCapacity = bCap;
		}
		
		static void setPermutationLength(int permLength)
		{
			permutationLength = permLength;
		}
		
		int permutationIndex(int index)
		{
			return permutation[index];
		}
		
		void setPermutationIndex(int index, int value)
		{
			permutation[index] = value;
		}
		
		int getNumberOfBins()
		{
			return numberOfBins;
		}
		
		vector<int>* getBinCapacities()
		{
			return &binCapacities;
		}
		
		vector<int>* getNumItemsPerBin()
		{
			return &numItemsPerBin;
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

int Chromosome::binCapacity;
int Chromosome::permutationLength;

class GA
{
	private:
		vector<Chromosome> population;
		int populationSize, numberOfGenerations, tournamentSize;
		double crossoverProbability, mutationProbability;
		int chromosomeLength;
		int binCap;
		int nBins;
	
	public:
		GA(int pSize, int ng, int tSize, double cP, double mP, int n, int b, vector<int> items)
		{
			populationSize = pSize;
			numberOfGenerations = ng;
			tournamentSize = tSize;
			crossoverProbability = cP;
			mutationProbability = mP;
			chromosomeLength = n;
			binCap = b;
			
			for (int counter = 0; counter < populationSize; counter++)
			{
				Chromosome::setPermutationLength(n);
				Chromosome::setBinCapacity(b);
				population.push_back(Chromosome(items));
			}
		}
		
		int getNBins()
		{
			return nBins;
		}
		
		static bool compare(Chromosome a, Chromosome b)
		{
			return (a.fitnessFunction() < b.fitnessFunction());
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
			fitness = best->fitnessFunction();
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
				double fit = tournament[counter]->fitnessFunction();
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
					double fit = tournament[counter]->fitnessFunction();
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
			double fitnessP1 = population[parents[0]].fitnessFunction();
			double fitnessP2 = population[parents[1]].fitnessFunction();
			double fitnessOffspring = offspring->fitnessFunction();
			
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
			double fitnessP1 = (*result).fitnessFunction();
			double fitnessOffspring = offspring->fitnessFunction();		
		
			if (fitnessP1 > fitnessOffspring)
			{
				population.erase(result);
				population.push_back(*offspring);
			}
		}
};

enum Type
{
	INT,
	DOUBLE,
	VECTOR_INT,
	ARITHMETHIC,
	SUMMATION
};

class TreeNode
{
	public:
		vector<void*> pointers;
		vector<Type> types;
		vector<TreeNode*> children;
		TreeNode()
		{
			
		}
		
		void printTreeNode()
		{
			cout << "\n";
			int size = pointers.size();
			for (int counter = 0; counter < size; counter++)
			{
				if (pointers[counter] != 0)
				{
					cout << "Not_Null - ";
				}
				cout << types[counter] << " - ";
				cout << "Child " << counter + 1 << ": ";
				if (children[counter] != 0)
				{
					children[counter]->printTreeNode();
				}
				cout << "\n";
			}
			cout << "\n";
		}
};

class GP
{
	public:
		
		TreeNode* root;
		int maxDepth;
	
		// Terminal Set
		vector<int> binCapacities;
		vector<int> numItemsPerBin;
		int binCapacity;
		int permutationLength;
	
		vector<void*> functionSet;
		vector<Type> functionTypes;
		vector<void*> terminalSet;
		vector<Type> terminalTypes;
	
		// Function Set
		static double plus (double arg1, double arg2)
		{
			return arg1 + arg2;
		}
		
		static double minus (double arg1, double arg2)
		{
			return arg1 - arg2;
		}
		
		static double multiply (double arg1, double arg2)
		{
			return arg1 * arg2;
		}
		
		static double divide (double arg1, double arg2)
		{
			return arg1 / arg2;
		}
		
		static double power (double arg1, double arg2)
		{
			double ans = 1;
			
			for (int count = 0; count < arg2; count++)
			{
				ans *= arg1;
			}
			
			return ans;
		}
		
		static double summation(vector<double>* values)
		{
			return std::accumulate(values->begin(), values->end(), 0.0);
		}
		
		GP(int maxD)
		{
			maxDepth = maxD;
			initialPopulationGeneration();
		}
		
		void initialPopulationGeneration()
		{
			srand(5);
			root = new TreeNode;
			functionSet.push_back((void*)&plus);
			functionTypes.push_back(ARITHMETHIC);
			functionSet.push_back((void*)&minus);
			functionTypes.push_back(ARITHMETHIC);
			functionSet.push_back((void*)&multiply);
			functionTypes.push_back(ARITHMETHIC);
			functionSet.push_back((void*)&divide);
			functionTypes.push_back(ARITHMETHIC);
			functionSet.push_back((void*)&power);
			functionTypes.push_back(ARITHMETHIC);
			functionSet.push_back((void*)&summation);
			functionTypes.push_back(SUMMATION);

			terminalSet.push_back((void*)&binCapacities);
			terminalTypes.push_back(VECTOR_INT);
			terminalSet.push_back((void*)&numItemsPerBin);
			terminalTypes.push_back(VECTOR_INT);
			terminalSet.push_back((void*)&binCapacity);
			terminalTypes.push_back(INT);
			terminalSet.push_back((void*)&permutationLength);
			terminalTypes.push_back(INT);
			
			root = gen_rnd_expr(root, maxDepth);
			root->printTreeNode();
		}
		
		TreeNode* gen_rnd_expr(TreeNode* node, int depth)
		{
			if (depth == 0)
			{
				int randomIndex = rand() % terminalSet.size();
				node->pointers.push_back(terminalSet[randomIndex]);
				node->types.push_back(terminalTypes[randomIndex]);
			}
			else
			{
				int randomIndex = rand() % functionSet.size();
				node->pointers.push_back(functionSet[randomIndex]);
				Type t = functionTypes[randomIndex];
				node->types.push_back(t);
				if (t == ARITHMETHIC)
				{
					for (int count = 0; count < 2; count++)
					{
						TreeNode* newNode = new TreeNode();

						node->children.push_back(gen_rnd_expr(newNode, depth - 1));
					}
				}
				else
				{
					int randomIndex = rand() % 20;
					for (int count = 0; count < randomIndex; count++)
					{
						TreeNode* newNode = new TreeNode();
						node->children.push_back(gen_rnd_expr(newNode, depth - 1));
					}
				}
			}
			return node;
		}
		
		void evolve ()
		{
			
		}
		
		double fitnessFunc ()
		{
			return 0.0;
		}
		
		void selection  ()
		{
			
		}

		void crossover  ()
		{
			
		}

		void mutation  ()
		{
			
		}

		void insertion  ()
		{
			
		}
		
};

int main ()
{
	/*
	ifstream inFile;
	string fileName;
	cin >> fileName;
	inFile.open(fileName.c_str());

	if (!inFile) 
	{
		cout << "Unable to open file.";
		exit(1);   // call system to stop
	}
	
	int counter = 0;
	int numberOfItems, binCapacity;
	vector<int> items;
	int x;
	while (inFile >> x) 
	{
		if (counter == 0)
		{
			numberOfItems = x;
		}
		else if (counter == 1)
		{
			binCapacity = x;
		}
		else
		{
			items.push_back(x);
		}
		counter++;
	}
	inFile.close();
	

	GA ga(500, 1000, 3, 1, 0.10, numberOfItems, binCapacity, items); 
	double fitness = ga.evolve();
	int numBins = ga.getNBins();
	*/
	
	GP (3);
	
	return 0;
}
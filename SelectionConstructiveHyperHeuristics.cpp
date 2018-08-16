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
#include <chrono>

using namespace std;
using namespace std::chrono;

// 1-Dimensional Bin-Packing Problem
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
		int totalBins;
		// Initialised by GA class.
		// 3 Training and Test Examples - One Per Level Of Difficulty
		int heuristicsLength;
		static vector<InputData>* input;
	
		/* Characteristics of the Problem:
		* Must be determined during the Packing - while applying different heuristics
		* Stored in these data structures:
		* Length Corresponds to number of bins for particular training/test instance
		*/
		vector<int> binCapacities;
		vector<int> numItemsPerBin;
		int sumItems;
		int sumBins;
		
		Chromosome()
		{
			sumItems = 0;
			sumBins = 0;
			totalBins = 0;
			heuristicsLength = (rand() % 15) + 1;
			// Random Initialisation of Fixed Length Non-repeating Heuristic Combinations
			for (int count = 0; count < heuristicsLength; count++)
			{
				heuristics.push_back(rand() % 4);
			}
			
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
			tabuListLength = 0;
			tempHeuristicsLength = 0;
			objectiveValueTabu = INT_MAX;
			totalBinsTabu = INT_MAX;
			copyLength = 0;
		}

		Chromosome(int heu)
        {
            sumItems = 0;
            sumBins = 0;
            totalBins = 0;
            heuristicsLength = 1;

            heuristics.push_back(heu);

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
            tabuListLength = 0;
            tempHeuristicsLength = 0;
            objectiveValueTabu = INT_MAX;
            totalBinsTabu = INT_MAX;
            copyLength = 0;
        }
		
		// Used When Creating Offspring
		void nonInitialisation(vector<int> h)
		{
			sumItems = 0;
			sumBins = 0;
			totalBins = 0;
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
			sumBins = 0;
			sumItems = 0;
			totalBins = 0;
			
			// Traverse Circular through the heuristic combination till items are exhausted.
			for (int count = 0; count < (*input)[num].numberOfItems; count++)
			{
				/* Pass the training/test instance number and the number of the current item 
				* being packed
				*/
				sumItems += (*input)[num].items[count];
				switch(heuristics[nextHeuristic])
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

		void checkPack(int num)
		{
			int calc = 0;
			for (int count = 0; count < numberOfBins[num]; count++)
			{
				calc += (*input)[num].binCapacity - binCapacities[count];
			}

			if (calc != sumItems)
			{
				cout << "PROBLEM" << "\n";
			}
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
		//	checkPack(num);
		//	cout << "";
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
		//	checkPack(num);
		//	cout << "";
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
		//	checkPack(num);
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
		//	checkPack(num);
		//	cout << "";
		}
		
		void fitnessFunction(int num)
		{
			double sum = 0.0;
			for (int count = 0; count < numberOfBins[num]; count++)
			{
				double calc = (*input)[num].binCapacity - binCapacities[count];
				sumBins += calc;
				calc /= (*input)[num].binCapacity;
				calc *= calc;
				
				sum += calc;
			}

			if (sumBins != sumItems)
			{
				cout << "WHAT\n";
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
				totalBins += numberOfBins[count];
			}
			
			return fit/3.0;
		}
		
		void printChromosome(bool train)
		{
			int start = 0;
			int end = 3;
			
			if (train != true)
			{
				start = 3;
				end = 6;
			}
			
			cout << "\nHeuristic Combination: ";
			for (int count = 0; count < heuristicsLength; count++)
			{
				cout << heuristics[count];
			}
			cout << " -- " << fitnesses[start] << ", " << fitnesses[start + 1] << ", " << fitnesses[start + 2];
			cout << "--";
			for (int k = start; k < end; k++)
			{
				cout << numberOfBins[k] << "--";
//				for (int count = 0; count < numberOfBins[k]; count++)
//				{
//					cout << binCapacities[count] << " ";
//				}
			}
			
			cout << "\n";
		}

		static void setInput(vector<InputData>* in)
		{
			input = in;
		}
		
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
		vector< vector<int> > tabuList;
		int tabuListLength;
		vector<int> tempHeuristics;
		int tempHeuristicsLength;
		vector<int> copy;
		int copyLength;
		static int tabuIterations;
		double objectiveValueTabu;
		int totalBinsTabu;
			
		template<typename T>
		bool isEqual(std::vector<T> const &v1, std::vector<T> const &v2)
		{
			return (v1.size() == v2.size() &&
			std::equal(v1.begin(), v1.end(), v2.begin()));
		}

		void tabu_search()
		{
			objectiveValueTabu = fitnessOverInputs(true);
			totalBinsTabu = totalBins;
			for (int count = 0; count < tabuIterations; count++)
			{
				moveOperator();
				bool tabu = false;
				for (int iterator = 0; iterator < tabuListLength; iterator++)
				{
					if (isEqual(tabuList[iterator], tempHeuristics))
					{
						tabu = true;
						break;
					}
				}
				
				if (tabu == false)
				{
					for (int w = 0; w < heuristicsLength; w++)
					{
						copy.push_back(heuristics[w]);
						copyLength++;
					}
					heuristics.clear();
					heuristicsLength = 0;
					
					for (int w = 0; w < tempHeuristicsLength; w++)
					{
						heuristics.push_back(tempHeuristics[w]);
						heuristicsLength++;
					}
					
					for (int g = 0; g < 3; g++)
					{
						applyHeuristics(g);
					}
					
					double ov = fitnessOverInputs(true);
					
					if (totalBins < totalBinsTabu || (totalBinsTabu == totalBins && ov < objectiveValueTabu))
					{
						objectiveValueTabu = ov;
						copyLength = 0;
						copy.clear();
					}
					else
					{
						heuristicsLength = 0;
						heuristics.clear();
						
						for (int w = 0; w < copyLength; w++)
						{
							heuristics.push_back(copy[w]);
							heuristicsLength++;
						}
						
						copyLength = 0;
						copy.clear();
					}
					
					vector<int>* addTabu = new vector<int>;
					for (int w = 0; w < tempHeuristicsLength; w++)
					{
						addTabu->push_back(tempHeuristics[w]);
					}
					
					tabuList.push_back(*addTabu);
					tempHeuristics.clear();
					tempHeuristicsLength = 0;
				}
			}	
		}
		
		void moveOperator()
		{
			int rand1 = rand() % 4;
			int rand2, rand3;
			switch(rand1)
			{
				case 0:
					TOO_LARGE:
					if (heuristicsLength == 1)
					{
						goto STOP_SEG_FAULT_;
					}
					rand2 = rand() % heuristicsLength;
					for (int w = 0; w < heuristicsLength; w++)
					{
						tempHeuristics.push_back(heuristics[w]);
					}
					tempHeuristics.erase(tempHeuristics.begin() + rand2);
					tempHeuristicsLength = heuristicsLength - 1;
					break;
				case 1:
					STOP_SEG_FAULT_:
					rand2 = rand() % 4;
					if (heuristicsLength > 8)
					{
						goto TOO_LARGE;
					}
					for (int w = 0; w < heuristicsLength; w++)
					{
						tempHeuristics.push_back(heuristics[w]);
					}
					tempHeuristics.push_back(rand2);
					tempHeuristicsLength = heuristicsLength;
					tempHeuristicsLength += 1;
					for (int w = 0; w < 3; w++)
					{
						rand2 = rand() % 4;
						tempHeuristics.push_back(rand2);
						tempHeuristicsLength += 1;
					}	
					break;
				case 2:
					rand2 = rand() % heuristicsLength;
					rand3 =  rand() % 4;
					if (heuristicsLength == 1)
					{
						goto STOP_SEG_FAULT_;
					}
					for (int w = 0; w < heuristicsLength; w++)
					{
						tempHeuristics.push_back(heuristics[w]);
					}
					while (rand3 == tempHeuristics[rand2])
					{
						rand3 =  rand() % 4;
					}
					tempHeuristics[rand2] = rand3;
					tempHeuristicsLength = heuristicsLength;
					break;
				case 3:
					rand2 = rand() % heuristicsLength;
					rand3 =  rand() % heuristicsLength;
					if (heuristicsLength == 1)
					{
						goto STOP_SEG_FAULT_;
					}
					while (rand3 == rand2)
					{
						rand3 =  rand() % heuristicsLength;
					}
					for (int w = 0; w < heuristicsLength; w++)
					{
						tempHeuristics.push_back(heuristics[w]);
					}
					tempHeuristics[rand2] = tempHeuristics[rand3];
					tempHeuristicsLength = heuristicsLength;
					break;
			}
		}
		
		static void setTabuIterations(int it)
		{
			tabuIterations = it;
		}
	
};

vector<InputData>* Chromosome::input;
int Chromosome::tabuIterations;

class GA
{
	private:
		vector<Chromosome> population;
		int populationSize, numberOfGenerations, tournamentSize;
		double crossoverProbability, mutationProbability;
		// Objective Value For Each Training and Test Instance
		vector<int>* nBins;
		// Best Chromosome Fitness Value
		Chromosome* best;
	
	public:
		GA(int pSize, int ng, int tSize, double cP, double mP, vector<InputData>* input)
		{
			populationSize = pSize;
			numberOfGenerations = ng;
			tournamentSize = tSize;
			crossoverProbability = cP;
			mutationProbability = mP;
			
			for (int counter = 0; counter < populationSize; counter++)
			{
				Chromosome::setInput(input);
				population.push_back(Chromosome());
				// population[counter].printChromosome(true);
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
			return ((a.totalBins < b.totalBins) || (a.totalBins == b.totalBins && a.fitnessOverInputs(true) < b.fitnessOverInputs(true)));
		}
		
		double evolve()
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

//			cout << "\nPopulation:";
//			for (int z = 0; z < populationSize; z++)
//			{
//				population[z].printChromosome(true);
//			}

			//cout << "\n\nBest Individual Training:";
			//best->printChromosome(true);
			
			return fitness; 
		}
		
		void selectParents(int p [2])
		{
			Chromosome* tournament[tournamentSize];
			vector<int> values;
			int valueLength = 0;
			
			for (int counter = 0; counter < tournamentSize; counter++)
			{
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
			int random = rand() % population[p[0]].heuristicsLength;
			int random2 = rand() % population[p[1]].heuristicsLength;

			vector<int> temp1;
			vector<int> temp2;

			for (int it = 0; it <= random; it++)
			{
				temp1.push_back(population[p[0]].heuristics[it]);
			}

			for (int it = 0; it <= random2; it++)
			{
				temp2.push_back(population[p[1]].heuristics[it]);
			}

			for (int it = random2 + 1; it < population[p[1]].heuristicsLength; it++)
			{
				temp1.push_back(population[p[1]].heuristics[it]);
			}

			for (int it = random + 1; it < population[p[0]].heuristicsLength; it++)
			{
				temp2.push_back(population[p[0]].heuristics[it]);
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
			//if ((rand() % 10000) >= (mutationProbability * 10000))
			//{
			//	return;
			//}
			int random1 = rand() % offspring->heuristicsLength;
			int random2 = rand() % 4;
			
			offspring->heuristics[random1] = random2;
			
			for (int count = 0; count < 3; count++)
			{
				offspring->applyHeuristics(count);
			}
		}
		
		void insert(int parents[2], Chromosome* offspring)
		{

			std::vector<Chromosome>::iterator result = std::max_element(population.begin(), population.end(), compare);

			double fitnessOffspring = offspring->fitnessOverInputs(true);
			double fitnessWorst = result->fitnessOverInputs(true);
			
			if (result->totalBins > offspring->totalBins || (result->totalBins == offspring->totalBins && fitnessWorst >= fitnessOffspring))
			{
					population.erase(result);
					population.push_back(*offspring);
			}
		}
};

int check_error_bits(ifstream* f) {
	int stop = 0;
	if (f->eof()) {
		perror("stream eofbit. error state");
		// EOF after std::getline() is not the criterion to stop processing
		// data: In case there is data between the last delimiter and EOF,
		// getline() extracts it and sets the eofbit.
		stop = 0;
	}
	if (f->fail()) {
		perror("stream failbit (or badbit). error state");
		stop = 1;
	}
	if (f->bad()) {
		perror("stream badbit. error state");
		stop = 1;
	}
	return stop;
}

int main ()
{
	ifstream inFile;
	ifstream names;
	string aName = "/home/vignesh/CLionProjects/SelectionConstructiveHyperHeuristics_Assignment1_COS790/selection-constructive-hyper-heuristics/FileNames.txt";
	names.open(aName.c_str());
	string temp = "";
	int u = 0;
	string fileName[6];
	while (1)
	{
		std::getline(names, temp);

		if (check_error_bits(&names))
		{
			break;
		}

		fileName[u] = temp;
		u++;
	}
	names.close();
	int counter = 0;
	vector<InputData> input;
	int x;

	for (int start = 0; start < 6; start++)
	{
		inFile.open(("/home/vignesh/CLionProjects/SelectionConstructiveHyperHeuristics_Assignment1_COS790/selection-constructive-hyper-heuristics/" + fileName[start]).c_str());
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

//	// GENETIC ALGORITHM
//	// Adjust Parameters
//	// Add More Parameters
//	int minTrainingTestBinSum = INT_MAX;
//	vector<double> bs;
//	vector<double> fits;
//	double avg = 0.0, standardDeviation = 0.0;
//	vector<double> fitnesses;
//	long avgDuration = 0;
//	for (int run = 0; run < 30; run++)
//	{
//		srand ((unsigned int) run + 1);
//		GA ga(100, 1000, 3, 1, 0.80, &input);
//		int nb = 0;
//
//		// Get starting timepoint
//		auto start = high_resolution_clock::now();
//		double fitnessTraining = ga.evolve();
//		auto stop = high_resolution_clock::now();
//		// Get ending timepoint
//
//		// Get duration. Substart timepoints to
//		// get durarion. To cast it to proper unit
//		// use duration cast method
//		auto duration = duration_cast<microseconds>(stop - start);
//
//		avgDuration += duration.count();
//
//		//cout << "Fitness Training: " << fitnessTraining << "\n";
//		vector<int> *numBinsTraining = ga.getNBins();
//		//cout << "Number of Bins Training: ";
//		for (int q = 0; q < 3; q++)
//		{
//		//	cout << (*numBinsTraining)[q] << " ";
//			nb += (*numBinsTraining)[q];
//		}
//		//cout << "\n";
//		Chromosome *trainingWinner = ga.getBest();
//
//		trainingWinner->applyHeuristics(3);
//		trainingWinner->applyHeuristics(4);
//		trainingWinner->applyHeuristics(5);
//		//cout << "\nBest Individual Testing: ";
//		cout << "ITERATION" << run;
//		trainingWinner->printChromosome(true);
//		trainingWinner->printChromosome(false);
//		cout << "____________________________________________\n";
//		double fitnessTesting = trainingWinner->fitnessOverInputs(false);
//		//cout << "Fitness Testing: " << fitnessTesting << "\n";
//		vector<int> *numBinsTesting = ga.getNBins();
//		//cout << "Number of Bins Testing: ";
//		for (int q = 3; q < 6; q++)
//		{
//		//	cout << (*numBinsTesting)[q] << " ";
//			nb += (*numBinsTesting)[q];
//		}
//
//		if (nb < minTrainingTestBinSum)
//		{
//			minTrainingTestBinSum = nb;
//			bs.clear();
//			for (int q = 0; q < 3; q++)
//			{
//				bs.push_back((*numBinsTraining)[q]);
//			}
//
//			for (int q = 3; q < 6; q++)
//			{
//				bs.push_back((*numBinsTesting)[q]);
//			}
//
//			fits.clear();
//			fits.push_back(fitnessTraining);
//			fits.push_back(fitnessTesting);
//			cout << "\n**********************************";
//			trainingWinner->printChromosome(true);
//			trainingWinner->printChromosome(false);
//			cout << "**********************************\n";
//		}
//		avg += fitnessTesting + fitnessTraining;
//		fitnesses.push_back(fitnessTesting);
//		fitnesses.push_back(fitnessTraining);
//		//cout << "\n\n\n";
//	}
//	cout << "\n";
//	avg /= 60;
//
//	for(int i = 0; i < 60; ++i)
//		standardDeviation += pow(fitnesses[i] - avg, 2);
//
//	for (int traverse = 0; traverse < 6; traverse++)
//	{
//		cout << bs[traverse] << " ";
//	}
//	cout << "\n";
//	for (int traverse = 0; traverse < 2; traverse++)
//	{
//		cout << fits[traverse] << " ";
//	}
//	cout << "\n";
//	cout << avg;
//	cout << "\n";
//	cout << sqrt(standardDeviation/60);
//	cout << "\n";
//	cout << avgDuration/30.0 * pow(10, -6);
//	cout << "\n";


//	// TABU LOCAL SEARCH
//	int minTrainingTestBinSum = INT_MAX;
//	vector<double> bs;
//	vector<double> fits;
//	double avg = 0.0, standardDeviation = 0.0;
//	vector<double> fitnesses;
//	long avgDuration = 0;
//	for (int run = 0; run < 30; run++)
//	{
//		srand ((unsigned int) run + 2);
//		Chromosome::setTabuIterations(10000);
//		Chromosome::setInput(&input);
//		Chromosome tabu;
//		int nb = 0;
//		// Get starting timepoint
//		auto start = high_resolution_clock::now();
//		tabu.tabu_search();
//		auto stop = high_resolution_clock::now();
//		// Get ending timepoint
//
//		// Get duration. Substart timepoints to
//		// get durarion. To cast it to proper unit
//		// use duration cast method
//		auto duration = duration_cast<microseconds>(stop - start);
//
//		avgDuration += duration.count();
//
//		double fitnessTraining = tabu.objectiveValueTabu;
//		//cout << "\n\n";
//		//cout << "Fitness Training: " << fitnessTraining << "\n";
//		vector<int> numBinsTraining = tabu.numberOfBins;
//		//cout << "Number of Bins Training: ";
//		for (int q = 0; q < 3; q++)
//		{
//		//		cout << (numBinsTraining)[q] << " ";
//			nb += (numBinsTraining)[q];
//		}
//		//cout << "\n";
//
//		tabu.applyHeuristics(3);
//		tabu.applyHeuristics(4);
//		tabu.applyHeuristics(5);
//
//		cout << "ITERATION" << run;
//		tabu.printChromosome(true);
//		tabu.printChromosome(false);
//		cout << "____________________________________________\n";
//
//		double fitnessTesting = tabu.fitnessOverInputs(false);
//		//cout << "Fitness Testing: " << fitnessTesting << "\n";
//		vector<int> numBinsTesting = tabu.numberOfBins;
//		//cout << "Number of Bins Testing: ";
//		for (int q = 3; q < 6; q++)
//		{
//		//		cout << (numBinsTesting)[q] << " ";
//			nb += (numBinsTesting)[q];
//		}
//		//cout << "\n\n\n";
//
//		if (nb < minTrainingTestBinSum)
//		{
//			minTrainingTestBinSum = nb;
//			bs.clear();
//			for (int q = 0; q < 3; q++)
//			{
//				bs.push_back((numBinsTraining)[q]);
//			}
//
//			for (int q = 3; q < 6; q++)
//			{
//				bs.push_back((numBinsTesting)[q]);
//			}
//
//			fits.clear();
//			fits.push_back(fitnessTraining);
//			fits.push_back(fitnessTesting);
//			cout << "\n**********************************";
//			tabu.printChromosome(true);
//			tabu.printChromosome(false);
//			cout << "**********************************\n";
//		}
//		avg += fitnessTesting + fitnessTraining;
//		fitnesses.push_back(fitnessTesting);
//		fitnesses.push_back(fitnessTraining);
//		//cout << "\n\n\n";
//	}
//	cout << "\n";
//	avg /= 60;
//
//	for(int i = 0; i < 60; ++i)
//		standardDeviation += pow(fitnesses[i] - avg, 2);
//
//	for (int traverse = 0; traverse < 6; traverse++)
//	{
//		cout << bs[traverse] << " ";
//	}
//	cout << "\n";
//	for (int traverse = 0; traverse < 2; traverse++)
//	{
//		cout << fits[traverse] << " ";
//	}
//	cout << "\n";
//	cout << avg;
//	cout << "\n";
//	cout << sqrt(standardDeviation/60);
//	cout << "\n";
//	cout << avgDuration/30.0 * pow(10, -6);
//	cout << "\n";

    srand (time(NULL));
	Chromosome::setInput(&input);
	Chromosome low_level_heuristic(3);
    low_level_heuristic.printChromosome(true);
    low_level_heuristic.applyHeuristics(3);
    low_level_heuristic.applyHeuristics(4);
    low_level_heuristic.applyHeuristics(5);
    low_level_heuristic.printChromosome(false);

    return 0;
}
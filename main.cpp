#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <random>
#include <boost/sort/sort.hpp>
#include <algorithm>    // std::sort
#include <thread>

#include "linearRegModel.h"

using namespace std;

bool compareModels (linearRegModel * model1, linearRegModel * model2);
void swapNodes(linearRegModelNode * node1, linearRegModelNode * node2);
void recombine(linearRegModel * parent1, linearRegModel * parent2, linearRegModel * child1, linearRegModel * child2, double probNodeSelection);
string getBest(vector<linearRegModel *> &modelArr);

int main()
{
    //Algorithm parameters
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    int nRuns = 30;
    int generations = 30;
    int population = 1000;
    int maxDepth = 5;
    int maxNodes = 50;
    int minDepth = 3;
    int tournamentSize = 2; //As a number of individuals
    int nRecombinations = 3;
    double terminalRate = 0.2;
    double constantRate = 0.2;
    double islandPercentSize = 0.01;
    double probNodeSelection = 0.01; //Probability of a node being selected for recombination
    double branchMutationRate = 0.8; //Number of branch mutations per generation per island as a percentage of the population per island
    double nodeMutationRate = 0.2;
    string filePath = "d1.csv";
    string outputFileName = "run_results.txt";


    //Program variable declarations and prologue
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    long double * rawData; //size allocated in file ingest
    int nRows; //size determined in file ingest
    int nCols; //size determined in file ingest
    int nIslands = (int)round(1/islandPercentSize);
    int islandSize = (int)round(population * islandPercentSize);

    string str; //Temp string. Used in file ingest

    vector<linearRegModel *> modelArr;
    vector<linearRegModel *> currentModelArr;
    vector<linearRegModel *> tournament;
    vector<thread> threadVector;

    ofstream results(outputFileName);

    time_t start, end; //Used to time the algorithm execution time

    srand(time(nullptr)); //rand must be seeded to create random models

    std::knuth_b rand_engine;  //random engine (knuth_b can be replace with another one)
    std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0); //Function to generate random double between 0.0 and 1.0


    //File ingest
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //open file
    ifstream fileData;
    fileData.open(filePath);

    //count number of rows to define array size (In case of variable file length)
    nRows = 0;
    while(getline(fileData, str))
    {
        nRows++;
    }
    //return to beginning of file
    fileData.clear();
    fileData.seekg(0, ios::beg);

    //count number of columns
    getline(fileData, str);
    nCols = 1;
    nCols += count(str.begin(), str.end(), ',');

    //return to beginning of file
    fileData.clear();
    fileData.seekg(0, ios::beg);

    //variable declarations needed for file ingest
    rawData = new long double[nRows*nCols];

    //file ingest (row by row)
    for(int i=0; i<nRows; i++)
    {
        for(int j=0; j<(nCols-1); j++)
        {
            getline(fileData, str, ',');
            *rawData = strtod(str.c_str(), nullptr);
            ++rawData;
        }
        getline(fileData, str);
        *rawData = strtod(str.c_str(), nullptr);
        ++rawData;
    }
    //return to beginning of rawData array
    rawData -= (nRows*nCols);

    //close file
    fileData.close();


    results << "ALGORITHM PARAMETERS" << endl;
    results << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    results << "Input File: " << filePath << endl;
    results << "nRuns: " << nRuns << endl;
    results << "generations: " << generations << endl;
    results << "population: " << population << endl;
    results << "maxDepth: " << maxDepth << endl;
    results << "maxNodes: " << maxNodes << endl;
    results << "minDepth: " << minDepth << endl;
    results << "tournamentSize: " << tournamentSize << endl;
    results << "nRecombinaitons: " << nRecombinations << endl;
    results << "terminalRate: " << terminalRate << endl;
    results << "constantRate: " << constantRate << endl;
    results << "islandPercentSize: " << islandPercentSize << endl;
    results << "probNodeSelection: " << probNodeSelection << endl;
    results << "branchMutationRate: " << branchMutationRate << endl;
    results << "nodeMutationRate: " << nodeMutationRate << endl;
    results << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    results << endl;
    results << "ALGORITHM RUNS" << endl;
    results << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    for(int run=0; run<nRuns; run++)
    {
        modelArr.reserve(population);
        tournament.reserve(tournamentSize);
        time(&start);
        //Random tree creation
        //------------------------------------------------------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < population; i++)
        {
            linearRegModel *tempModel = new linearRegModel(nCols, maxDepth, maxNodes, minDepth, terminalRate, constantRate);
            modelArr.emplace_back(tempModel);
            modelArr[i]->calculateRootMeanSquaredError(rawData, nRows);
        }

        //Run Algorithm
        //------------------------------------------------------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < generations; i++) {
            currentModelArr.clear();
            currentModelArr.reserve(round(population + (population * (branchMutationRate * 2))));

            //copy modelArr into currentModelArr
            for(int j = 0; j < population; j++) {
                linearRegModel *tempModel = new linearRegModel;
                *tempModel = *modelArr[j];
                currentModelArr.emplace_back(tempModel);
            }

            //Node Mutation
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //mutate X% of the models but keep best model in each island
            for (int j = 0; j < nIslands; j++)
            {
                //eliteism
                linearRegModel *tempModel = new linearRegModel;
                *tempModel = **min_element(modelArr.begin() + (j * islandSize),
                                           modelArr.begin() + ((j + 1) * islandSize),
                                           compareModels);
                for (int k = 0; k < islandSize; k++) {
                    //Mutate
                    if (nodeMutationRate >= uniform_zero_to_one(rand_engine))
                    {
                        //Generate new thread
                        static thread t1([&]()
                        {
                            currentModelArr[(j * islandSize) + k]->mutateNode(currentModelArr[(j * islandSize) + k]->getRandomNode(probNodeSelection));
                            currentModelArr[(j * islandSize) + k]->calculateRootMeanSquaredError(rawData, nRows);
                        });
                        threadVector.push_back(move(t1));
                    }
                }
                //worst is replaced by best
                **max_element(modelArr.begin() + (j * islandSize),
                              modelArr.begin() + ((j + 1) * islandSize),
                              compareModels)
                        = *tempModel;
            }
            //join all threads
            for(int t=0; t<threadVector.size(); t++)
            {
                if(threadVector[t].joinable())
                    threadVector.at(t).join();
            }
            threadVector.clear();

            //Branch Mutation
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //Select parents from island and create children from them
            //loop through each island
            for (int j = 0; j < nIslands; j++)
            {
                //loop through number of branch mutations
                for (int k = 0; k < (int) round(islandSize * branchMutationRate); k++)
                {
                    //Generate new thread
                    //static thread t([&]()
                    //{
                        //select parents
                        linearRegModel *parent1 = currentModelArr[j * islandSize + rand() % islandSize];
                        linearRegModel *parent2 = new linearRegModel(nCols, maxDepth, maxNodes, minDepth, terminalRate,
                                                                     constantRate); //Diversity injection

                        //create children
                        linearRegModel *child1 = new linearRegModel;
                        linearRegModel *child2 = new linearRegModel;

                        //recombine
                        recombine(parent1, parent2, child1, child2, probNodeSelection);

                        //recalculate standard deviation
                        child1->calculateRootMeanSquaredError(rawData, nRows);
                        child2->calculateRootMeanSquaredError(rawData, nRows);

                        //append children to array
                        currentModelArr.emplace_back(child1);
                        currentModelArr.emplace_back(child2);
                    //});
                    //threadVector.push_back(move(t));
                }
                //join all threads
                //for(int t=0; t<threadVector.size(); t++)
                //{
                //    if(threadVector[t].joinable())
                //        threadVector.at(t).join();
                //}
                //threadVector.clear();
            }

            //Rebuild Islands
            int currentIslandSize = (int) round(islandSize + (islandSize * branchMutationRate * 2));
            //loop through islands
            for (int j = (nIslands - 1); j >= 0; j--) {
                //loop through children per island
                for (int k = 0; k < (int) round(islandSize * branchMutationRate * 2); k++) {
                    //insert at beginning of island
                    linearRegModel *tempModel = new linearRegModel;
                    *tempModel = *currentModelArr.back();
                    currentModelArr.pop_back();
                    currentModelArr.insert(currentModelArr.begin() + (j * islandSize), tempModel);
                }
            }

            //Shuffle Islands
            for (int j = 0; j < nIslands; j++) {
                shuffle(currentModelArr.begin() + (int) round(j * currentIslandSize),
                        currentModelArr.begin() + (int) round((j + 1) * currentIslandSize),
                        rand_engine);
            }

            //Survivor Selection
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //Make modelArr ready to receive survivors
            modelArr.clear();
            modelArr.reserve(population);

            //Use tournament selection by island
            for (int j = 0; j < nIslands; j++) {
                //eliteism
                linearRegModel *tempModel = new linearRegModel;
                *tempModel = **min_element(currentModelArr.begin() + (j * currentIslandSize),
                                           currentModelArr.begin() + ((j + 1) * currentIslandSize),
                                           compareModels);
                for (int k = 0; k < (islandSize - 1); k++) {
                    for (int t = 0; t < tournamentSize; t++) {
                        tournament[t] = currentModelArr[j * currentIslandSize + rand() % currentIslandSize];
                    }
                    linearRegModel *tempModel = new linearRegModel;
                    *tempModel = **min_element(tournament.begin(), tournament.end(), compareModels);
                    modelArr.emplace_back(tempModel);
                }
                modelArr.emplace_back(tempModel);
            }

            //Inter-Island recombination
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //Shift all elements over by one
            //modelArr.insert(modelArr.begin(), modelArr.back());
            //modelArr.pop_back();

            currentModelArr.clear();
            currentModelArr.reserve(round(population + (population * (branchMutationRate * 2))));

            //Sort models in islands to find best parents for recombination
            for (int j = 0; j < nIslands; j++) {
                boost::sort::spinsort(modelArr.begin() + (j * islandSize),
                                      modelArr.begin() + ((j + 1) * (islandSize)),
                                      compareModels);
            }

            //copy modelArr into currentModelArr
            for (int j = 0; j < population; j++) {
                linearRegModel *tempModel = new linearRegModel;
                *tempModel = *modelArr[j];
                currentModelArr.emplace_back(tempModel);
            }

            //Recombine with best parents from two separate islands
            for (int j = 0; j < nIslands; j++) {
                for (int k = 0; k < nRecombinations; k++) {
                    //select best parents (first element from each island)
                    linearRegModel *parent1 = currentModelArr[j * islandSize + k];
                    linearRegModel *parent2 = currentModelArr[((j + 1) * islandSize + k) % population];

                    //create children
                    linearRegModel *child1 = new linearRegModel;
                    linearRegModel *child2 = new linearRegModel;

                    //recombine
                    recombine(parent1, parent2, child1, child2, probNodeSelection);

                    //recalculate standard deviation
                    child1->calculateRootMeanSquaredError(rawData, nRows);
                    child2->calculateRootMeanSquaredError(rawData, nRows);

                    //append children to array
                    currentModelArr.emplace_back(child1);
                    currentModelArr.emplace_back(child2);
                }
            }

            //rebuild islands
            for (int j = (nIslands - 1); j >= 0; j--) {
                for (int k = 0; k < (nRecombinations * 2); k++) {
                    linearRegModel *tempModel = new linearRegModel;

                    //insert at beginning of island
                    *tempModel = *currentModelArr.back();
                    currentModelArr.insert(currentModelArr.begin() + ((j * islandSize) + (j * nRecombinations * 2)),
                                           tempModel);
                    currentModelArr.pop_back();
                }
            }

            //shuffle islands
            for (int j = 0; j < nIslands; j++) {
                shuffle(currentModelArr.begin() + ((j * islandSize) + (j * nRecombinations * 2)),
                        currentModelArr.begin() + (((j + 1) * islandSize) + ((j + 1) * nRecombinations * 2)),
                        rand_engine);
            }

            //Survivor Selection
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //made modelArr ready to receive survivors
            modelArr.clear();
            modelArr.reserve(population);

            //Use tournament selection by island
            currentIslandSize = islandSize + (nRecombinations * 2);
            for (int j = 0; j < nIslands; j++) {
                //eliteism
                linearRegModel *tempModel = new linearRegModel;
                *tempModel = **min_element(currentModelArr.begin() + (j * currentIslandSize),
                                           currentModelArr.begin() + ((j + 1) * currentIslandSize),
                                           compareModels);
                for (int k = 0; k < (islandSize - 1); k++) {
                    for (int t = 0; t < tournamentSize; t++) {
                        tournament[t] = currentModelArr[j * currentIslandSize + rand() % currentIslandSize];
                    }
                    linearRegModel *tempModel = new linearRegModel;
                    *tempModel = **min_element(tournament.begin(), tournament.end(), compareModels);
                    modelArr.emplace_back(tempModel);
                }
                modelArr.emplace_back(tempModel);
            }

            //Determine fitness
            //------------------------------------------------------------------------------------------------------------------------------------------------------------
            //print best standard deviation

            cout << "Generation " << i << ":" << endl;
            //cout << getBest(modelArr);
            //cout << endl;

            if(i == generations-1)
            {
                time(&end);
                cout << "writing to file" << endl;
                results << "Run: " << run << endl
                        << getBest(modelArr) << endl
                        << "Execution time: " << double(end - start) << " seconds" << endl
                        << endl << endl;
            }
        }
        modelArr.clear();
        tournament.clear();
    }
    results.close();
    return 0;
}

bool compareModels (linearRegModel * model1, linearRegModel * model2)
{
    return (model1->getRootMeanSquaredError() < model2->getRootMeanSquaredError());
}

void recombine(linearRegModel * parent1, linearRegModel * parent2, linearRegModel * child1, linearRegModel * child2, double probNodeSelection)
{
    child1->copyModel(parent1);
    child2->copyModel(parent2);

    swapNodes(child1->getRandomNode(probNodeSelection), child2->getRandomNode(probNodeSelection));
}

void swapNodes(linearRegModelNode * node1, linearRegModelNode * node2)
{

    swap(node1->data, node2->data);
    swap(node1->terminal, node2->terminal);

    //Swap leftChild
    if(node1->leftChild == nullptr)
    {
        if(node2->leftChild != nullptr)
        {
            node1->leftChild = node2->leftChild;
            node2->leftChild = nullptr;
        }
        //else both are nullptr so nothing needs to be done
    }
    else if(node2->leftChild == nullptr)
    {
        node2->leftChild = node1->leftChild;
        node1->leftChild = nullptr;
    }
    else
    {
        swap(node1->leftChild, node2->leftChild);
    }

    //Swap rightChild
    if(node1->rightChild == nullptr)
    {
        if(node2->rightChild != nullptr)
        {
            node1->rightChild = node2->rightChild;
            node2->rightChild = nullptr;
        }
        //else both are nullptr so nothing needs to be done
    }
    else if(node2->leftChild == nullptr)
    {
        node2->rightChild = node1->rightChild;
        node1->rightChild = nullptr;
    }
    else
    {
        swap(node1->rightChild, node2->rightChild);
    }
}

string getBest(vector<linearRegModel *> &modelArr)
{
    stringstream ss;
    linearRegModel * tempModel;
    tempModel = *min_element(modelArr.begin(), modelArr.end(), compareModels);
    ss << "RMSE: " << tempModel->getRootMeanSquaredError() << endl
       << "Model: " << tempModel->getModel();
    return ss.str();
}
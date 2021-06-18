#ifndef LINEARREGMODEL_H
#define LINEARREGMODEL_H

//#include <cstdlib>
#include <iostream>
#include <random>

struct linearRegModelNode
{
    bool terminal;
    std::string data;
    linearRegModelNode * leftChild;
    linearRegModelNode * rightChild;
};

class linearRegModel
{
    public:
        linearRegModel() = default;
        linearRegModel(int nCols_param, int maxDepth_param, int maxNodes_param, int minDepth_param, double terminalRate_param, double constantRate_param);
        virtual ~linearRegModel();

        int getMinDepth(){return minDepth;}
        int getMaxDepth(){return maxDepth;}
        int getMaxNodes(){return maxNodes;}
        int getSize(){return nCurrentNodes;}
        int getnCols(){return nCols;}
        double getTerminalRate(){return terminalRate;}
        double getConstantRate(){return constantRate;}
        double getRootMeanSquaredError(){return rootMeanSquaredError;}
        linearRegModelNode * getRoot(){return root;}

        void showModel(){showModelNode(root);}
        std::string getModel(){return getModelNode(root);}

        void calculateRootMeanSquaredError(long double * rawData, int nRows);
        void mutateRandomNode(){parseNodesMutate(root, rand()%nCurrentNodes);}
        void copyModel(linearRegModel * parent);

        linearRegModelNode * getRandomNode(double probRecombination);
        void mutateNode(linearRegModelNode * currentNode);

    private:
        linearRegModelNode * root;
        int currentDepth;
        int nCurrentNodes;
        int minDepth;
        int maxDepth;
        int maxNodes;
        int nCols;
        bool error;
        double terminalRate;
        double constantRate;
        double rootMeanSquaredError;

        std::knuth_b rand_engine;
        //std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);

        linearRegModelNode * makeTerminal(linearRegModelNode * currentNode);
        void showModelNode(linearRegModelNode * currentNode);
        linearRegModelNode * makeModel();
        long double calculateExpected(long double * rawData, linearRegModelNode * currentNode);
        int parseNodesMutate(linearRegModelNode * currentNode, int nodesToVisit);
        linearRegModelNode * copyModelNode(linearRegModelNode * parent);
        linearRegModelNode * parseRandomNode(linearRegModelNode * currentNode, double probRecombination);
        std::string getModelNode(linearRegModelNode * currentNode);

        //int parseNodesRecombine(linearRegModelNode * parent1, linearRegModelNode * parent2, linearRegModelNode * child1, linearRegModelNode * child2, int nodesToVisit);

};

#endif // LINEARREGMODEL_H

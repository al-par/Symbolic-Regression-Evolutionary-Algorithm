//Class Header
#include "linearRegModel.h"

//Included libraries
#include <math.h>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <boost/algorithm/string.hpp>

//Function macros
#define LEFT_VAL calculateExpected(rawData, currentNode->leftChild)
#define RIGHT_VAL calculateExpected(rawData, currentNode->rightChild)


//Parameterized constructor
linearRegModel::linearRegModel(int nCols_param, int maxDepth_param, int maxNodes_param, int minDepth_param, double terminalRate_param, double constantRate_param)
{
    bool fail = false;
    if(nCols_param <= 1)
    {
        std::cout << "Data must contain at least one dependent and one independent variable. Number of columns must be greater than or equal to 2." << std::endl;
        fail = true;
    }
    if((maxDepth_param < 1) || (maxDepth_param < minDepth_param))
    {
        std::cout << "Max depth [" << maxDepth_param << "] is too small. Must be greater than 0 and greater than min depth." << std::endl;
        fail = true;
    }
    if(maxNodes_param < 1)
    {
        std::cout << "Max nodes [" << maxDepth_param << "] is too small. Must be greater than 0." << std::endl;
        fail = true;
    }
    if(terminalRate_param <= 0)
    {
        std::cout << "Warning: terminal rate is [" << terminalRate_param << "]. Model will generate until max nodes [" << maxNodes_param << "] is reached." << std::endl;
    }
    if(constantRate_param <= 0)
    {
        std::cout << "Warning: constant rate is [" << constantRate_param << "]. Model will not generate any constants." << std::endl;
    }
    if(fail)
    {
        //std::cout << "Returning blank linearRegModel." << std::endl;
        std::cout << "Destructing linearRegModel." << std::endl;
        this->~linearRegModel();
    }
    else
    {
        nCols = nCols_param;
        minDepth = minDepth_param;
        maxDepth = maxDepth_param;
        maxNodes = maxNodes_param;
        terminalRate = terminalRate_param;
        constantRate = constantRate_param;

        error = false;
        currentDepth = 0;
        nCurrentNodes = 0;

        root = makeModel();
    }
}

linearRegModel::~linearRegModel()
{
    //dtor
}

linearRegModelNode * linearRegModel::makeModel()
{
    ++currentDepth;
    ++nCurrentNodes;
    linearRegModelNode * currentNode = new linearRegModelNode;

    //node is terminal
    if((((rand()%100) <= (terminalRate*100)) && (currentDepth >= minDepth)) || (currentDepth >= maxDepth) || (nCurrentNodes >= maxNodes))
    {
        currentNode = makeTerminal(currentNode);
    }
    else //node is an operator
    {
        currentNode->terminal = false;
        //select operator
        int x = rand()%9;
        switch(x)
        {
            case 0: currentNode->data = "+"; break;
            case 1: currentNode->data = "-"; break;
            case 2: currentNode->data = "/"; break;
            case 3: currentNode->data = "*"; break;
            case 4: currentNode->data = "sin"; break;
            case 5: currentNode->data = "cos"; break;
            case 6: currentNode->data = "log10"; break;
            case 7: currentNode->data = "log"; break;
            case 8: currentNode->data = "exp"; break;
            //case 9: currentNode->data = "pow"; break;
        }

        //create children
        //operation does not require two operands
        //if((x >= 4) && (x < 9))
        if(x >= 4)
        {
            currentNode->leftChild = makeModel();
            currentNode->rightChild = nullptr;
        }
        else //operation requires two operands
        {
            currentNode->leftChild = makeModel();
            currentNode->rightChild = makeModel();
        }
    }
    --currentDepth;
    return currentNode;
}

linearRegModelNode * linearRegModel::makeTerminal(linearRegModelNode * currentNode)
{
    currentNode->terminal = true;
    currentNode->leftChild = nullptr;
    currentNode->rightChild = nullptr;

    //terminal is a constant
    mutateNode(currentNode);
    return currentNode;
}

void linearRegModel::showModelNode(linearRegModelNode * currentNode)
{
    //Node is a variable or constant
    if(currentNode->terminal)
    {
        std::cout << " " << currentNode->data << " ";
    }
    else //Node is operator
    {
        //Operation takes two operands
        if(currentNode->rightChild)
        {
            std::cout << " (";
            showModelNode(currentNode->leftChild);
            std::cout << ") " << currentNode->data << " (";
            showModelNode(currentNode->rightChild);
        }
        else //Operation takes one operand
        {
            std::cout << " " << currentNode->data << " ( ";
            showModelNode(currentNode->leftChild);
        }
        std::cout << ") ";
    }
}

void linearRegModel::calculateRootMeanSquaredError(long double * rawData, int nRows)
{
    long double * rawDataStart = rawData;

    rootMeanSquaredError = 0;
    for(int i=0; i<nRows; i++)
    {
        //(expected - actual)^2
        rootMeanSquaredError += pow((calculateExpected(rawData, root) - rawData[nCols-1]), 2);
        rawData += nCols;
        if(error) //restart calculation
        {
            i = 0;
            rootMeanSquaredError = 0;
            rawData = rawDataStart;
            error = false;
        }
    }
    rootMeanSquaredError = sqrt(rootMeanSquaredError/nRows); //Root mean squared error

    //If standard deviation is undefined or infinity, then this model is useless. Try a new model
    if(isinf(rootMeanSquaredError) ||  isnan(rootMeanSquaredError))
    {
        root = makeModel();
        calculateRootMeanSquaredError(rawDataStart, nRows);
    }
}

long double linearRegModel::calculateExpected(long double * rawData, linearRegModelNode * currentNode)
{
    std::string str = currentNode->data;
    //current node is a terminal
    if(currentNode->terminal)
    {
        //terminal is a variable
        if(str.find('c') != str.npos)
        {
            boost::erase_all(str, "c");
            return rawData[std::stoi(str)];
        }
        else //Terminal is a constant
        {
            return strtod(str.c_str(), nullptr);
        }
    }
    else //current node is an operator
    {
        while(true)
        {
            if (str == "+")          {return (LEFT_VAL + RIGHT_VAL);}
            else if (str == "-")     {return (LEFT_VAL - RIGHT_VAL);}
            else if (str == "/")     {long double x = RIGHT_VAL; if(x != 0){return (LEFT_VAL/x);}} //Check that we don't divide by zero
            else if (str == "*")     {return (LEFT_VAL * RIGHT_VAL);}
            else if (str == "sin")   {return sin(LEFT_VAL);}
            else if (str == "cos")   {return cos(LEFT_VAL);}
            else if (str == "log10") {long double x = LEFT_VAL; if(x > 0){return log10(x);}} //Check that we don't perform an illegal operaion
            else if (str == "log")   {long double x = LEFT_VAL; if(x > 0){return log(x);}} //Check that we don't perform an illegal operaion
            else if (str == "exp")   {return exp(LEFT_VAL);}
            //else if (str == "pow")   {return pow(LEFT_VAL, RIGHT_VAL);}

            //No legal operation was possible. Mutate the node and restart the calculation
            mutateNode(currentNode);
            str = currentNode->data;
            error = true;
        }
    }
}

void linearRegModel::mutateNode(linearRegModelNode * currentNode)
{
    //Node is a terminal
    if(currentNode->leftChild == nullptr)
    {
        if((rand()%100) <= (constantRate*100))
        {
            currentNode->data = std::to_string((double)((rand()%100)) / ((double)((rand()%100)+1)));
        }
        else //terminal is a variable
        {
            currentNode->data = "c";
            currentNode->data.append(std::to_string(rand()%(nCols-1)));
        }
    }
    else if(currentNode->rightChild == nullptr) //Node is operator with only one operand
    {
        switch(rand()%5)
        {
            case 0: currentNode->data = "sin"; break;
            case 1: currentNode->data = "cos"; break;
            case 3: currentNode->data = "log10"; break;
            case 4: currentNode->data = "log"; break;
            case 5: currentNode->data = "exp"; break;
        }
    }
    else //Node is operator with two operands
    {
        switch(rand()%4)
        {
            case 0: currentNode->data = "+"; break;
            case 1: currentNode->data = "-"; break;
            case 2: currentNode->data = "*"; break;
            case 3: currentNode->data = "/"; break;
            //case 4: currentNode->data = "pow"; break;
        }
    }
}

int linearRegModel::parseNodesMutate(linearRegModelNode * currentNode, int nodesToVisit)
{
    nodesToVisit--;
    if(nodesToVisit < 0)
    {
        return nodesToVisit;
    }
    else if(nodesToVisit == 0) //Mutate this node
    {
        //std::cout << std::endl << "NODE MUTATED" << std::endl << "Before: ";
        //showModel();
        //std::cout << std::endl << "Mutation: " << currentNode->data << " --> ";
        mutateNode(currentNode);
        //std::cout << currentNode->data << std::endl;
        //std::cout << "After: ";
        //showModel();
        return nodesToVisit;
    }
    else if(currentNode->terminal) //Node is a terminal
    {
        return nodesToVisit;
    }
    else //Node is operator
    {
        nodesToVisit = parseNodesMutate(currentNode->leftChild, nodesToVisit);
        if(nodesToVisit < 0)
        {
            return nodesToVisit;
        }
        else if(currentNode->rightChild != nullptr)
        {
            return parseNodesMutate(currentNode->rightChild, nodesToVisit);
        }
        return nodesToVisit;
    }
}

linearRegModelNode * linearRegModel::getRandomNode(double probSelection)
{
    linearRegModelNode * tempNode = nullptr;
    while(tempNode == nullptr)
    {
        tempNode = parseRandomNode(root, probSelection);
    }
    return tempNode;
}

linearRegModelNode * linearRegModel::parseRandomNode(linearRegModelNode * currentNode, double probSelection)
{
    std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);
    linearRegModelNode * tempNode;

    //Current node is selected
    if(probSelection >= uniform_zero_to_one(rand_engine))
    {
        return currentNode;
    }
    else if(currentNode->leftChild != nullptr) //Node is not a terminal and not is selected
    {
        tempNode = parseRandomNode(currentNode->leftChild, probSelection);
        if(tempNode != nullptr)
        {
            return tempNode; //Node was selected
        }
        else if(currentNode->rightChild != nullptr) //Node has right child and is not selected
        {
            tempNode = parseRandomNode(currentNode->rightChild, probSelection);
            if(tempNode != nullptr)
            {
                return tempNode;
            }
        }
    }
    return nullptr; //Node is not selected
}

void linearRegModel::copyModel(linearRegModel * parent)
{
    currentDepth;
    nCurrentNodes;
    minDepth = parent->getMinDepth();
    maxDepth = parent->getMaxDepth();
    maxNodes = parent->getMaxNodes();
    nCols = parent->getnCols();
    terminalRate = parent->getTerminalRate();
    constantRate = parent->getConstantRate();
    root = copyModelNode(parent->getRoot());
}

linearRegModelNode * linearRegModel::copyModelNode(linearRegModelNode * parent)
{
    linearRegModelNode * currentChildNode = new linearRegModelNode;

    if(parent->terminal)
    {
        currentChildNode->leftChild = nullptr;
        currentChildNode->rightChild = nullptr;
    }
    else
    {
        currentChildNode->leftChild = copyModelNode(parent->leftChild);
        if(parent->rightChild)
        {
            currentChildNode->rightChild = copyModelNode(parent->rightChild);
        }
        else
        {
            currentChildNode->rightChild = nullptr;
        }
    }
    currentChildNode->terminal = parent->terminal;
    currentChildNode->data = parent->data;
    return currentChildNode;
}

std::string linearRegModel::getModelNode(linearRegModelNode * currentNode)
{
    std::stringstream ss;
    //Node is a variable or constant
    if(currentNode->terminal)
    {
        ss  << currentNode->data;
    }
    else //Node is operator
    {
        //Operation takes two operands
        if(currentNode->rightChild)
        {
            ss << "(";
            ss << getModelNode(currentNode->leftChild);
            ss << ")" << currentNode->data << "(";
            ss << getModelNode(currentNode->rightChild);
        }
        else //Operation takes one operand
        {
            ss  << currentNode->data << "(";
            ss << getModelNode(currentNode->leftChild);
        }
        ss << ")";
    }
    return ss.str();
}

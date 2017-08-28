#ifndef _FRAMEWORK
#define _FRAMEWORK

#include <vector>
#include <map>
#include <string>

#include "Matrix.H"
#include "NCALearner.H"
#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Graph.H"
#include "GraphLearner.H"

#include <map>
#include <vector>
using namespace std;

class Framework
{
public:
	Framework();
	~Framework();
	int init(int argc, char** argv);
	int start();
	int resume();
	bool shouldResume();
private:
	bool lShouldResume;
	bool lResumeMERLIN;
	int  resIter;//resume iteration
	bool doNextStep(int itrNum, EdgeList*& newNet);
	int readTFA(int itrNum);
	VariableManager* regMngr;
	VariableManager* tgtMngr;
	VariableManager* allGenesMngr;
	EvidenceManager* evidMngr;
	EvidenceManager* allEvidMngr;
	Graph* priorNet;
	EdgeList* pnet;
	double lambda; //LASSO penalty
	int CV;//CV
	int numNcaRep;
	double sparsity;
	double beta;
	int numThread;
	char outdir[1025];

	int printHelp(char* name);

	GraphLearner* getGL(int itrNum, int subNum, map<string,int>& geneModuleID, EdgeList* initNet, bool sl);

	static void* runOneGL(void* ptr);
	
	EdgeList* inferNetwork(int itrNum, EdgeList* initNet);
	int inferNCA(int itrNum, EdgeList* initNet);

	bool isConverged(int itrNum, EdgeList* oldNet, EdgeList* newNet);

	map<string,int> moduleID;
	int readClusters(char* aFName,map<string,int>& geneModuleID);
	int writeClusters(char* aFName,map<string,int>& geneModuleID);
};

#endif //_FRAMEWORK

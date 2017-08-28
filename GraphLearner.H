#ifndef _GRAPHLEARNER
#define _GRAPHLEARNER

#include <map>
#include "Potential.H"
#include "FactorGraph.H"
#include "CommonTypes.H"
#include "EdgeList.H"
using namespace std;


//This is gamma map for each variable within each datapoint
typedef map<int,INTDBLMAP*> DATA_GAMMA_MAP;
typedef map<int,INTDBLMAP*>::iterator DATA_GAMMA_MAP_ITER;

class FactorManager;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class GraphLearner
{
	public:
	GraphLearner();
	~GraphLearner();
	int setInputFName(const char*);
	int setMaxFactorSize(int);
	int setMaxFactorSize_Approx(int);
	int setPenalty(double);
	int setBeta1(double aval);
	int setBeta_Motif(double bval);
	int setLambda(double bval);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setRestrictedList(VariableManager* tempMgr);
	int setRestrictedList(vector<string>* names);
	int setPreRandomizeSplit();
	int setTargetList(const char*);
	int setNoPrune(bool);
	int init();
	int initPartitions(int);
	int initEdgePriorMeta();
	int initEdgePriorMeta(map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePrior);
	int doCrossValidation(int);
	int start(int);
	int start_gradualMBIncrease(int);
	int start_gradualMBIncrease_RankRegulators(int);
	int getTopRegs(map<int,int>&);
	int setGlobalEvidenceManager(EvidenceManager*);
	int setHoldOutEvManager(EvidenceManager*);
	int setVariableManager(VariableManager*);
	int setOutputDirName(const char*);
	int setClusteringThreshold(double);
	int setSpecificFold(int);
	int showModelParameters(int);
	int getPredictionError_Holdout(int);
	int getPredictionError_CrossValid(int);
	int setRandom(bool);
	int readModuleMembership(const char*);
	int setDefaultModuleMembership();

	//int setMotifGraph(const char*);
	//int setChIPGraph(const char*);
	//int setBeta_ChIP(double bval);
	int setPriorGraph_All(const char*);
	int setDefaultModuleMembership(VariableManager* tempVarMngr);

	int setDefaultModuleMembership(vector<string>* names);
	int setDefaultModuleMembership(map<string,int>& moduleID);
	int setPriorGraph_All(const char* netname, EdgeList* pnet, double b);
	EdgeList* getLearnedGraph();
	map<string,int>* getLearnedModule();

	int setInitNet(EdgeList* i_n);

	int setShouldLoad(bool f);
	int setTopOutputDirName(const char* dirPath);
	private:
	
	//map<int,INTDBLMAP*> edgePriors_Motif;
	//map<int,INTDBLMAP*> edgePriors_ChIP;
	//map<string,map<string,double>*> priorGraph_Motif;
	//map<string,map<string,double>*> priorGraph_ChIP;
	//double beta_chip;
	//int initEdgePriorMeta_Motif();
	//int initEdgePriorMeta_ChIP();
	int initEdgePriorMeta_All();
	map <string,double> betamap;
	map <string,map<int,INTDBLMAP*>*> edgepriormap;
	map <string,map<string,map<string,double>*>*> priorgraphmap;
	map <string,EdgeList*> priorpointers;
		
	double getInitPLLScore();
	double getPLLScore();
	double getInitPrior();
	double getPriorChange();
	double getEdgePrior(int tfID, int targetID);
	double getPLLScore_Datapoint(int datapoint);
	int clearFoldSpecData();
	int initEdgeSet(bool);
	int initVariableComponents();
	int initCondsetMap_Nopool();
	INTINTMAP* getConditionSet(int cind);

	int collectMoves(int);
	int collectMoves(int,int);
	int sortMoves();
	int makeMoves();
	int attemptMove(MetaMove*,map<int,INTINTMAP*>&);
	int collectRemoveMoves(int,int);
	int makeRemoveMoves();
	int attemptRemoveMove(MetaMove*,map<int,INTINTMAP*>&);
	int getPLLScoreRemoveOne(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot);
	int getNewPLLScore(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement,Potential**);
	double getNewPLLScore_Condition(int cId, int csetId, int vId, Potential* newPot,int);
	double getNewPLLScore_Condition(int csetId, int vId, Potential* newPot);
	double getNewPLLScore_Condition(int csetId, int vId, int uID, Potential* newPot);
	double getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uID, Potential* newPot);
	double getNewPLLScore_Condition_Noprior(int csetId, int vId, Potential* newPot);
	double getMBPrior(int targetID);

	int genCondSetKey(INTINTMAP& condSet, string& aKey);

	int cleanUp();
	int dumpAllGraphs(int,int,int);
	int getEdgeVars(string&, char*, char*);
	bool checkMBSize(INTINTMAP*,int,int,int);
	bool checkMBSize(int,int,int,int);
	int populateGraphsFromFile(const char* aFName);
	int populateGraphsFromInput();
	
	int writeIterCnt(int currFold, int iter);
	int readIterCnt(const char* aFName);
	int doTar();
	int readAddRemEdges(const char* aFName, map<string,int>& emap);
	int writeAddRemEdges(const char* aFName, map<string,int>& emap);
	map<string,int> resumeEdges;
	double loadInitPLLScore();
	int writePLLScore();
	int loadLastUpdate();
	int writeLastUpdate();



	int normalizeWeight(INTDBLMAP&);
	//Update the potentials to reflect the new covariance and mean estimates and also recompute the likelihood.
	int updatePotentials();

	int setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph);
	int redefineModules();
	int redefineModules_Global();
	double getModuleContrib(string& tgtName, string& tfName);
	double getModuleContribLogistic(string& tgtName, string& tfName);
	int initPhysicalDegree();
	int getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet);
	double getModuleWideScoreImprovement(int cid,INTINTMAP& conditionSet,Variable* u, Variable* v,map<string,int>* moduleMembers);
	
	map<int,FactorManager*> fgMgrSet;
	map<int,PotentialManager*> potMgrSet;
	map<int,EvidenceManager*> evMgrSet;
	map<int,FactorGraph*> fgGraphSet;

	EvidenceManager* globalEvMgr;
	EvidenceManager* holdoutEvMgr;

	map<int,string> outLocMap;

	VariableManager* varManager;
	int maxFactorSize;
	int maxFactorSizeApprox;
	int maxMBCntGreedy;
	char inputFName[1024];
	double penalty;
	double convThreshold;
	char restrictedFName[1024];
	map<string,INTINTMAP*> edgeConditionMap;
	FactorGraph* fg;
	vector<MetaMove*> moveSet;
	vector<MetaMove*> removeMoveSet;
	map<string,double> edgeValue;
	map<int,INTINTMAP*> condsetMap;
	map<string,int> condsetKeyIDMap;
	map<int,string> condsetIDKeyMap;
	map<string,int> restrictedVarList;
	map<int,INTDBLMAP*> currPLLMap;
	char trueGraphFName[1024];
	map<int,double> finalScores;
	map<int,INTINTMAP*> affectedVariables;
	bool preRandomizeSplit;
	gsl_rng* rnd;
	char outputDirName[1024];
	char topOutputDirName[1024];

	map<int,INTINTMAP*> datasetInd;
	map<string,double> edgePresenceProb;
	map<int,double> varNeighborhoodPrior;
	//These are the edges that got added
	map<string,int> addedEdges;
	map<string,int> removedEdges;
	map<string,int> edgeUpdates;
	map<string,int> edgeRemoveUpdates;
	double beta1;
	double beta_motif;
	double lambda;
	map<string,int> reqdTargetList;
	bool random;
	bool randomData;
	bool noprune;
	map<int,map<string,int>*> moduleGeneSet;
	map<string,int> geneModuleID;
	map<int,map<string,int>*> moduleIndegree;
	map<string,int> regulatorModuleOutdegree;
	map<int,int> idVidMap;
	int globalID;
	double clusterThreshold;
	map<string,int> variableStatus;
	int currFold;
	char foldoutDirName[1024];
	bool shouldLoad;
	int specificFold;
	EdgeList* initnet;
};
#endif //_GRAPHLEARNER

#ifndef _FACTOR_GRAPH
#define _FACTOR_GRAPH

#include <map>
#include <fstream>
#include "CommonTypes.H"
#include "EdgeList.H"
using namespace std;
class Variable;
class SlimFactor;
class VariableManager;
class FactorManager;

class FactorGraph
{
	public:
	FactorGraph();
	~FactorGraph();
	int setFactor(SlimFactor*);
	int getFactorCnt();
	SlimFactor* getFactorAt(int);
	map<int,SlimFactor*>& getAllFactors();
	int dumpFactors_ClusterFormat(const char*,int,VSET& variableSet);
	int dumpFactors_PairwiseFormat(const char*,int,VSET& variableSet);
	int dumpVarMB_ClusterFormat(const char*,int,VSET& variableSet,FactorManager*);
	int dumpVarMB_PairwiseFormat(const char*, int, VSET&,FactorManager*);
	int dumpVarMB_PairwiseFormat(ofstream&, VSET&);
	int dumpBias(ofstream&, VSET&,map<string,int>&);
	int readVarMB_PairwiseFormat(const char* aFName,VariableManager* vMgr,FactorManager* fMgr);
	int dumpCandidateVarMB_PairwiseFormat(const char*,int,VSET&,FactorManager*);
	bool isConsistent();
	int  getGraph(EdgeList* onet,VSET& variableSet);
	private:
	int showFactor_PW(SlimFactor*,ostream&,VSET& variableSet);
	map<int,SlimFactor*> factorSet;
};
#endif

#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "LatticeStructure.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "FactorManager.H"
#include "MetaMove.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "HyperGeomPval.H"
#include "Distance.H"
#include "GraphLearner.H"

GraphLearner::GraphLearner()
{
	shouldLoad = false;
	convThreshold=0.0001;
	restrictedFName[0]='\0';
	trueGraphFName[0]='\0';
	preRandomizeSplit=false;
	random=false;
	noprune=false;
	lambda=0;
	clusterThreshold=0.5;
	specificFold=-1;
	holdoutEvMgr=NULL;
	initnet = NULL;
	removedEdges.clear();
	addedEdges.clear();
	priorpointers.clear();
	priorgraphmap.clear();
	resumeEdges.clear();
}

GraphLearner::~GraphLearner()
{
	if(globalEvMgr != NULL)
	{
		delete globalEvMgr;
	}
	if (initnet != NULL)
	{
		delete initnet;
	}
	if (priorpointers.size()>0)
	{
		for (auto itr=priorpointers.begin();itr!=priorpointers.end();itr++)
		{
			EdgeList* e = itr->second;
			delete e;
		}
		priorpointers.clear();
	}
	for (auto itr=priorgraphmap.begin();itr!=priorgraphmap.end();itr++)
	{
		auto tmap = itr->second;
		for (auto titr=tmap->begin();titr!=tmap->end();titr++)
		{
			delete titr->second;
		}
		tmap->clear();
	}
	edgepriormap.clear();
	for (auto itr=edgepriormap.begin();itr!=edgepriormap.end();itr++)
	{
		auto tmap = itr->second;
		for (auto titr=tmap->begin();titr!=tmap->end();titr++)
		{
			delete titr->second;
		}
		tmap->clear();
	}
	edgepriormap.clear();
	clearFoldSpecData();
	for(map<int,FactorManager*>::iterator fIter=fgMgrSet.begin();fIter!=fgMgrSet.end();fIter++)
	{
		FactorManager* temp = fIter->second;
		delete temp;
	}
	varNeighborhoodPrior.clear();
	edgePresenceProb.clear();
}

int
GraphLearner::setShouldLoad(bool f)
{
	shouldLoad = f;
	return 0;
}

int
GraphLearner::setInitNet(EdgeList* i_n)
{
	initnet = i_n;
	return 0;
}

int
GraphLearner::setInputFName(const char* aFName)
{
	strcpy(inputFName,aFName);
	return 0;
}

int
GraphLearner::setMaxFactorSize(int aVal)
{
	maxFactorSize=aVal;
	return 0;
}

int
GraphLearner::setMaxFactorSize_Approx(int aVal)
{
	maxFactorSizeApprox=aVal;
	return 0;
}

int 
GraphLearner::setPenalty(double aVal)
{
	penalty=aVal;
	return 0;
}

int
GraphLearner::setBeta1(double aval)
{
	beta1=aval;
	return 0;
}

/*
int
GraphLearner::setBeta_ChIP(double aval)
{
	beta_chip=aval;
	return 0;
}

int
GraphLearner::setMotifGraph(const char* aFName)
{
	setPriorGraph(aFName,priorGraph_Motif);	
	return 0;
}

int
GraphLearner::setChIPGraph(const char* aFName)
{
	setPriorGraph(aFName,priorGraph_ChIP);
	return 0;
}

int
GraphLearner::initEdgePriorMeta_Motif()
{
	initEdgePriorMeta(priorGraph_Motif,edgePriors_Motif);
	return 0;
}

int
GraphLearner::initEdgePriorMeta_ChIP()
{	
	initEdgePriorMeta(priorGraph_ChIP,edgePriors_ChIP);
	return 0;
}

*/

int
GraphLearner::initEdgePriorMeta_All()
{
	for(map<string,map<string,map<string,double>*>*>::iterator gIter=priorgraphmap.begin();gIter!=priorgraphmap.end();gIter++)
	{
		map<string,map<string,double>*>* priorgraph = gIter->second;
		map<int,INTDBLMAP*>* edgeprior = new map<int,INTDBLMAP*>();
		edgepriormap[gIter->first] = edgeprior;
		initEdgePriorMeta(*priorgraph,*edgeprior);
	}
	return 0;
}

int
GraphLearner::setPriorGraph_All(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string gname;
		string fname;
		double gbeta;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				gname.append(tok);
			}
			else if(tokCnt==1)
			{
				fname.append(tok);
			}
			else if(tokCnt==2)
			{
				gbeta=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		betamap[gname] = gbeta;
		map<string,map<string,double>*>* priorGraph = new map<string,map<string,double>*>();
		setPriorGraph(fname.c_str(),*priorGraph);
		priorgraphmap[gname] = priorGraph;
	}
	inFile.close();
	return 0;
}

int
GraphLearner::setBeta_Motif(double aval)
{
	beta_motif=aval;
	return 0;
}

int
GraphLearner::setLambda(double l)
{
	lambda=l;
	return 0;
}

int 
GraphLearner::setConvergenceThreshold(double aVal)
{
	convThreshold=aVal;
	return 0;
}

int
GraphLearner::setRestrictedList(const char* aFName)
{
	strcpy(restrictedFName,aFName);
	ifstream inFile(restrictedFName);
	string buffer;
	while(inFile.good())
	{
		getline(inFile,buffer);
		if(buffer.length()<=0)
		{
			continue;
		}
		restrictedVarList[buffer]=0;
	}
	inFile.close();
	return 0;
}


int 
GraphLearner::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int
GraphLearner::setGlobalEvidenceManager(EvidenceManager* anEvMgr)
{
	globalEvMgr=anEvMgr;
	return 0;
}

int
GraphLearner::setHoldOutEvManager(EvidenceManager* aMgr)
{	
	holdoutEvMgr=aMgr;
	return 0;
}


int 
GraphLearner::setVariableManager(VariableManager* aPtr)
{
	varManager=aPtr;
	return 0;
}

int
GraphLearner::setTopOutputDirName(const char* dirPath)
{
	strcpy(topOutputDirName,dirPath);
	return 0;
}

int
GraphLearner::setOutputDirName(const char* dirPath)
{
	strcpy(outputDirName,dirPath);
	char mkdirCmd[1024];
	sprintf(mkdirCmd,"mkdir -p %s",outputDirName);
	system(mkdirCmd);
	return 0;
}

int 
GraphLearner::setClusteringThreshold(double aVal)
{
	clusterThreshold=aVal;
	return 0;
}
 

int
GraphLearner::setSpecificFold(int fid)
{
	specificFold=fid;
	return 0;
}

int 
GraphLearner::setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tfName;
		string tgtName;
		double edgeStrength;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgtName.append(tok);
			}
			else if(tokCnt==2)
			{
				edgeStrength=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,double>* tgtSet=NULL;
		if(priorGraph.find(tfName)==priorGraph.end())
		{
			tgtSet=new map<string,double>;
			priorGraph[tfName]=tgtSet;
		}
		else
		{
			tgtSet=priorGraph[tfName];
		}
		(*tgtSet)[tgtName]=edgeStrength;
	}
	inFile.close();
	return 0;
}


int 
GraphLearner::setTargetList(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tgtname(buffer);
		reqdTargetList[tgtname]=0;
	}
	inFile.close();
	return 0;
}

int
GraphLearner::setRandom(bool flag)
{
	random=flag;
	return 0;
}

int
GraphLearner::readModuleMembership(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		int moduleID;
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else if(tokCnt==1)
			{
				moduleID=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* geneSet=NULL;
		if(moduleGeneSet.find(moduleID)==moduleGeneSet.end())
		{
			geneSet=new map<string,int>;
			moduleGeneSet[moduleID]=geneSet;
		}
		else
		{
			geneSet=moduleGeneSet[moduleID];
		}
		(*geneSet)[geneName]=0;
		geneModuleID[geneName]=moduleID;
	}
	inFile.close();
	return 0;
}

int
GraphLearner::setDefaultModuleMembership()
{
	VSET& varSet=varManager->getVariableSet();
	int vCnt=varSet.size();
	int moduleCnt=(int) sqrt(vCnt/2);
	if(moduleCnt>30)
	{
		moduleCnt=30;
	}
	map<int,int> matIdvIdMap;
	int mID=0;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		matIdvIdMap[mID]=vIter->first;
		mID++;
	}
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//Randomly partition the variables into clusterassignments
	vector<int> randIndex;
	double step=1.0/(double)vCnt;
	map<int,int> usedInit;
	int maxind=0;
	for(int i=0;i<vCnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
			if(rind>maxind)
			{
				maxind=rind;
			}
		}
		//if(matIdvIdMap.find(rind)==matIdvIdMap.end())
		//{
		//	cout <<"Did not find " << rind << " matrixidmap " << endl;
		//}
		int dataind=matIdvIdMap[rind];
		usedInit[rind]=0;
		randIndex.push_back(dataind);
	}
	//For each partition estimate the mean and covariance
	int clusterSize=vCnt/moduleCnt;
	for(int e=0;e<moduleCnt;e++)
	{
		int startInd=e*clusterSize;
		int endInd=(e+1)*clusterSize;
		if(e==moduleCnt-1)
		{
			endInd=clusterSize;
		}
		map<string,int>* geneSet=NULL;
		geneSet=new map<string,int>;
		moduleGeneSet[e]=geneSet;
		for(int i=startInd;i<endInd;i++)
		{
			int dataId=randIndex[i];
			Variable* v=varSet[dataId];
			(*geneSet)[v->getName()]=0;
			geneModuleID[v->getName()]=e;
		}
	}
	randIndex.clear();
	matIdvIdMap.clear();
	usedInit.clear();
	return 0;
}

int 
GraphLearner::setNoPrune(bool pruneStatus)
{
	noprune=pruneStatus;
	return 0;
}

int
GraphLearner::initPartitions(int numberOfComponents)
{
	int randseed=0;
	globalEvMgr->partitionData(numberOfComponents,evMgrSet,randseed,datasetInd);
	char outputLoc[1024];
	char commandName[1024];

	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		EvidenceManager* evMgr=eIter->second;
		int datasetId=eIter->first;
		PotentialManager* potMgr=new PotentialManager;
		potMgr->setEvidenceManager(evMgr);
		evMgr->setVariableManager(varManager);
		//sprintf(outputLoc,"%s/model%d/",outputDirName,datasetId);
		sprintf(outputLoc,"%s",outputDirName);
		//sprintf(commandName,"mkdir %s",outputLoc);
		//system(commandName);
		string outputLocKey(outputLoc);
		outLocMap[datasetId]=outputLocKey;
		potMgr->setOutputDir(outputLoc);
		potMgrSet[datasetId]=potMgr;
	}
	cout <<"Created data for " << potMgrSet.size() << " partitions" << endl;
	return 0;
}

int
GraphLearner::initEdgePriorMeta(map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePriors)
{
	VSET& varSet=varManager->getVariableSet();
	for(map<string,int>::iterator rIter=restrictedVarList.begin();rIter!=restrictedVarList.end();rIter++)
	{
		int regId=varManager->getVarID(rIter->first.c_str());
		if(regId==-1)
		{
			continue;
		}
		if(graph.find(rIter->first)==graph.end())
		{
			continue;
		}
		int tfhit=0;
		map<string,double>* tgtSet=graph[rIter->first];
		for(map<string,double>::iterator vIter=tgtSet->begin();vIter!=tgtSet->end();vIter++)
		{
			INTDBLMAP* edgePriorGene=NULL;
			int tgtId=varManager->getVarID(vIter->first.c_str());
			if(tgtId==-1)
			{
				continue;
			}
			if(edgePriors.find(tgtId)==edgePriors.end())
			{
				edgePriorGene=new INTDBLMAP;
				edgePriors[tgtId]=edgePriorGene;
			}
			else
			{
				edgePriorGene=edgePriors[tgtId];
			}
			double ewt=fabs(vIter->second);
			if(edgePriorGene->find(regId)==edgePriorGene->end())
			{
				//(*edgePriorGene)[regId]=vIter->second;
				(*edgePriorGene)[regId]=ewt;
			}
			else
			{
				//(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+vIter->second;
				(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+ewt;
			}
			tfhit++;
		}
		//cout <<"TF: "<< rIter->first << " TFhits= " << tfhit << endl;
	}

	return 0;
}

int
GraphLearner::doCrossValidation(int foldCnt)
{
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	rnd=gsl_rng_alloc(gsl_rng_default);
	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		eIter->second->setFoldCnt(foldCnt);
		eIter->second->splitData(0);
	}
	//The first key is for the fold number
	//For each fold we have a trained model. For each trained model we have the likelihood on 
	//all the test sets, including the self test.
	int foldBegin=0;
	int foldEnd=foldCnt;
	if(specificFold>-1)
	{
		foldBegin=specificFold;
		foldEnd=specificFold+1;
	}
	for(int f=foldBegin;f<foldEnd;f++)
	{	
		for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
		{
			EvidenceManager* evMgr=eIter->second;
			evMgr->splitData(f);
			if(random)
			{	
				evMgr->randomizeEvidence(r);
			}
			PotentialManager* potMgr=potMgrSet[eIter->first];
			potMgr->reset();
			potMgr->resetCache();
			potMgr->setRandom(random);
			potMgr->init(f);
			if(fgMgrSet.find(eIter->first)!=fgMgrSet.end())
			{
				FactorManager* delMe=fgMgrSet[eIter->first];
				delete delMe;
			}
			FactorManager* fgMgr=new FactorManager;
			fgMgr->setPotentialManager(potMgr);
			fgMgr->setEvidenceManager(evMgr);
			fgMgr->setVariableManager(varManager);
			fgMgr->setOutputDir(outLocMap[eIter->first].c_str());
			fgMgr->setMaxFactorSize(maxFactorSize);
			fgMgr->setMaxFactorSize_Approx(maxFactorSizeApprox);
			fgMgr->setPenalty(penalty);
			fgMgrSet[eIter->first]=fgMgr;
			//if(strlen(restrictedFName)>0)
			//{
			//	fgMgr->readRestrictedVarlist(restrictedFName);
			//}
			fgMgr->setRestrictedList(restrictedVarList);
			fgMgr->allocateFactorSpace();

			if(fgMgr->readRandomInfo()==-1)
			{
				//cout <<"Did not find random mutual informations. Calling estimateRandomInfo " << endl;
				fgMgr->estimateRandomInfo_Approximate(SAMPLE_CNT);
				fgMgr->readRandomInfo();
			}
			if(fgMgr->readStructure(f)==-1)
			{
				fgMgr->learnStructure();
				fgMgr->showStructure_allK(f);
			}
			char outputDir[1024];
			sprintf(outputDir,"%s/fold%d",outLocMap[eIter->first].c_str(),f);
			char foldOutputDirCmd[1024];
			sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
			system(foldOutputDirCmd);
		}
		clearFoldSpecData();
		//start(f);
		start_gradualMBIncrease(f);
		//start_gradualMBIncrease_RankRegulators(f);
		int totalEvid=0;
		if(holdoutEvMgr!=NULL)
		{
			totalEvid=holdoutEvMgr->getNumberOfEvidences();
		}
		getPredictionError_CrossValid(f);
		if(totalEvid>0)
		{
			getPredictionError_Holdout(f);
		}
	}
	gsl_rng_free(r);

	/*
	char scoreFName[1024];
	sprintf(scoreFName,"%s/scoreFile.txt",outLocMap.begin()->second.c_str());
	ofstream sFile(scoreFName);
	for(map<int,double>::iterator dIter=finalScores.begin();dIter!=finalScores.end();dIter++)
	{
		sFile<< dIter->second << endl;
	}
	sFile.close();
	*/

	gsl_rng_free(rnd);
	return 0;
}

int
GraphLearner::start(int f)
{
	//Repeat until convergence
	//int currK=1;
	currFold=f;
	int maxMBSizeApprox=maxFactorSizeApprox-1;
	int currK=maxMBSizeApprox;
	rnd=gsl_rng_alloc(gsl_rng_default);
	//initEdgePriorMeta_Motif();
	//initEdgePriorMeta_ChIP();
	initEdgePriorMeta_All();
	initEdgeSet(false);
	initPhysicalDegree();
	int i=0;
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		idVidMap[i]=vIter->first;
		i++;
	}
	if(strlen(trueGraphFName)==0)
	{
		
		double currGlobalScore=getInitPLLScore();
		double initScore=getInitPrior();
		//currGlobalScore=currGlobalScore+initScore;
		int showid=0;
		int moduleiter=0;
		bool notConvergedTop=true;
		double scorePremodule=currGlobalScore;
		while(moduleiter<10 && notConvergedTop)
		{
			int subiter=0;
			while(subiter<varSet.size())
			{
				bool notConverged=true;
				int iter=0;
				int attemptedMoves=0;
				//while(notConverged && subiter<6000)
				while(notConverged && iter<200)
				{
					collectMoves(currK,subiter);
					if(moveSet.size()==0)
					{
						notConverged=false;
						continue;
					}
					sortMoves();
					makeMoves();
					double newScore=getPLLScore();
					//double newImprScore=getPriorChange();
					//newScore=newScore+newImprScore;
					double diff=newScore-currGlobalScore;
					if(diff<=convThreshold)
					{
						notConverged=false;
					}
					dumpAllGraphs(currK,f,showid);
					currGlobalScore=newScore;
					//cout <<"Current iter " << iter << " Score after beta-theta " << newScore << endl;
					for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
					{
						cIter->second->clear();
						delete cIter->second;
					}
					iter++;
					affectedVariables.clear();
					showid++;
					attemptedMoves++;
				}
				subiter++;
			}
			if(moduleiter>0)
			{
				if((currGlobalScore-scorePremodule)<convThreshold)
				{
					notConvergedTop=false;
				}
			}
			scorePremodule=currGlobalScore;
			//redefineModules();
			redefineModules_Global();
			moduleiter++;
		}
		cout <<"Final Score " << currGlobalScore << endl;
		finalScores[f]=currGlobalScore;
	}
	return 0;
}


int
GraphLearner::start_gradualMBIncrease(int f)
{
	int removemoves=0;
	int addmoves=0;
	//Repeat until convergence
	//int currK=1;
	currFold=f;
	sprintf(foldoutDirName,"%s/fold%d",outLocMap.begin()->second.c_str(),f);
	int maxMBSizeApprox=maxFactorSizeApprox-1;
	int currK=maxMBSizeApprox;
	rnd=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(rnd,rseed);
	//cout <<rseed << endl;
	//initEdgePriorMeta_Motif();
	//initEdgePriorMeta_ChIP();
	initEdgePriorMeta_All();
	initEdgeSet(false);
	initPhysicalDegree();
	int i=0;
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		idVidMap[i]=vIter->first;
		i++;
	}
	if (initnet != NULL)
	{
		//cout << "Populate from init!" << endl;
		if (!shouldLoad)
		{
			populateGraphsFromInput();
		}
		//makeInitMoves();
	}
	if(strlen(trueGraphFName)==0)
	{
		double currGlobalScore = 0;
		int titer=-1;
		if (shouldLoad)
		{
			cout << "Read modules..." << endl;
			string& dirname=outLocMap[evMgrSet.begin()->first];
			char mname[1024];
			sprintf(mname,"%s/fold%d/modules.txt",dirname.c_str(),f);
			readModuleMembership(mname);
			cout << "Read networks..." << endl;
			char fname[1024];
			sprintf(fname,"%s/prediction_k%d.txt",foldoutDirName,currK+1);
			populateGraphsFromFile(fname);
			char aoname[1024];
			sprintf(aoname,"%s/fold%d/added.txt",dirname.c_str(),f);
			readAddRemEdges(aoname,addedEdges);
			char roname[1024];
			sprintf(roname,"%s/fold%d/removed.txt",dirname.c_str(),f);
			readAddRemEdges(roname,removedEdges);
			cout << "Read iterations..." << endl;
			char iname[1024];
			sprintf(iname,"%s/fold%d/iter.txt",dirname.c_str(),f);
			titer = readIterCnt(iname);
			cout << "TITER:" << titer << endl;
			//for(map<string,int>::iterator mitr=geneModuleID.begin(); mitr!=geneModuleID.end();mitr++)
			//{
			//	variableStatus[mitr->first]=titer;
			//}
			currGlobalScore=loadInitPLLScore();
			loadLastUpdate();
		}
		else
		{
			currGlobalScore=getInitPLLScore();
		}

		double initScore=getInitPrior();
		//currGlobalScore=currGlobalScore+initScore;
		int showid=0;
		int moduleiter=0;
		bool notConvergedTop=true;
		vector<int> randOrder;
		while(moduleiter<1 && notConvergedTop)
		{
			int iter=titer+1;
			bool notConverged=true;
			while(notConverged && iter<50)
			{
				int attemptedMoves=0;
				int subiter=0;
				double scorePremodule=currGlobalScore;
				EvidenceManager* evMgr=evMgrSet.begin()->second;
				randOrder.clear();
				evMgr->populateRandIntegers(rnd,randOrder,varSet.size(),varSet.size());				
				struct timeval begintime;
				struct timeval endtime;
				struct timezone begintimezone;
				struct timezone endtimezone;
				gettimeofday(&begintime,&begintimezone);
				while(subiter<varSet.size())
				//while(notConverged && subiter<6000)
				{
					int rID=randOrder[subiter];
					if(idVidMap.find(rID)==idVidMap.end())
					{
						cout <<"Variable at  " << rID << " just not found " << endl;
						exit(0);
					}
					//int vID=idVidMap[rID];
					int vID=idVidMap[subiter];
					VSET_ITER vIter=varSet.find(vID);
					if(vIter==varSet.end())
					{
						subiter++;
						continue;
					}
					Variable* v=varSet[vID];
					int lastiter=0;
					if(variableStatus.find(v->getName())!=variableStatus.end())
					{
						lastiter=variableStatus[v->getName()];
						if((iter-lastiter)>=5)
						{
							//cout <<"Skipping " << v->getName() << endl;
							subiter++;
							continue;	
						}
					}
					else
					{
						variableStatus[v->getName()]=iter;
					}
					//collectMoves(currK,subiter);
					struct timeval begintime_v;
					struct timeval endtime_v;
					struct timezone begintimezone_v;
					struct timezone endtimezone_v;
					collectRemoveMoves(currK,vID);
					if(removeMoveSet.size()>0)
					{
						removemoves++;
						makeRemoveMoves();
					}
					collectMoves(currK,vID);
					if(moveSet.size()==0)
					{
						subiter++;
						continue;
					}
					else
					{
						addmoves++;
					}
					sortMoves();
					makeMoves();
					//double newScore=getPLLScore();
					//double newImprScore=getPriorChange();
					//newScore=newScore+newImprScore;
					//double diff=newScore-currGlobalScore;
					//if(diff<=convThreshold)
					//{
					//	notConverged=false;
					//}
					//dumpAllGraphs(currK,f,iter);
					//currGlobalScore=newScore;
					//cout <<"Current iter " << iter << " Score after beta-theta " << newScore << endl;
					for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
					{
						cIter->second->clear();
						delete cIter->second;
					}
					subiter++;
					affectedVariables.clear();
					showid++;
					attemptedMoves++;
					gettimeofday(&endtime_v,&endtimezone_v);
					//printf("Time elapsed for one var %uj secs %d microsec\n",(unsigned int)(endtime_v.tv_sec-begintime_v.tv_sec),(unsigned int)(endtime_v.tv_usec-begintime_v.tv_usec));
				}
				gettimeofday(&endtime,&endtimezone);
				double newScore=getPLLScore();
				currGlobalScore=newScore;
				//printf("Time elapsed for all vars %d mins %d secs %d microsec\n", (unsigned int)(endtimezone.tz_minuteswest-begintimezone.tz_minuteswest), (unsigned int)(endtime.tv_sec-begintime.tv_sec,endtime.tv_usec-begintime.tv_usec));
				if((currGlobalScore-scorePremodule)<=convThreshold)
				{
					notConverged=false;
				}
				else
				{
					redefineModules_Global();
				}

				scorePremodule=currGlobalScore;
				dumpAllGraphs(currK,f,iter);
				writeIterCnt(f,iter);
				writePLLScore();
				writeLastUpdate();
				string& dirname=outLocMap[evMgrSet.begin()->first];
				char aoname[1024];
				sprintf(aoname,"%s/fold%d/added.txt",dirname.c_str(),f);
				writeAddRemEdges(aoname,addedEdges);
				char roname[1024];
				sprintf(roname,"%s/fold%d/removed.txt",dirname.c_str(),f);
				writeAddRemEdges(roname,removedEdges);
				doTar();
				iter++;
			}
			moduleiter++;
		}
		cout <<"Final Score " << currGlobalScore << endl;
		finalScores[f]=currGlobalScore;
	}
	return 0;
}


int
GraphLearner::start_gradualMBIncrease_RankRegulators(int f)
{
	//Repeat until convergence
	//int currK=1;
	int maxMBSizeApprox=maxFactorSizeApprox-1;
	int currK=maxMBSizeApprox;
	rnd=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(rnd,rseed);
	//cout <<rseed << endl;
	//initEdgePriorMeta_Motif();
	//initEdgePriorMeta_ChIP();
	initEdgePriorMeta_All();
	initEdgeSet(false);
	initPhysicalDegree();
	int i=0;
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		idVidMap[i]=vIter->first;
		i++;
	}
	if(strlen(trueGraphFName)==0)
	{
		
		double currGlobalScore=getInitPLLScore();
		double initScore=getInitPrior();
		//currGlobalScore=currGlobalScore+initScore;
		int showid=0;
		int moduleiter=0;
		bool notConvergedTop=true;
		int k=1;
		while(moduleiter<1 && notConvergedTop)
		{
			int iter=0;
			bool notConverged=true;
			while(notConverged && iter<200)
			{
				int attemptedMoves=0;
				int subiter=0;
				double scorePremodule=currGlobalScore;
				bool notConvergedPerK=true;
				//while(subiter<varSet.size())
				while(notConvergedPerK)
				{
					collectMoves(k);
					if(moveSet.size()==0)
					{
						notConvergedPerK=false;	
						continue;
					}
					map<int,int> topRegs;
					getTopRegs(topRegs);
					map<int,MetaMove*> keepMoves;
					for(int i=0;i<moveSet.size();i++)
					{
						MetaMove* m=moveSet[i];
						if(topRegs.find(m->getSrcVertex())!=topRegs.end())
						{
							keepMoves[i]=m;
						}
						else
						{
							delete m;
						}
					}
					moveSet.clear();
					for(map<int,MetaMove*>::iterator mIter=keepMoves.begin();mIter!=keepMoves.end();mIter++)
					{
						moveSet.push_back(mIter->second);
					}
					keepMoves.clear();
					makeMoves();
					topRegs.clear();
					double newScore=getPLLScore();
					//double newImprScore=getPriorChange();
					//newScore=newScore+newImprScore;
					double diff=newScore-currGlobalScore;
					if(diff<=convThreshold)
					{
						notConvergedPerK=false;
					}
					dumpAllGraphs(currK,f,showid);
					currGlobalScore=newScore;
					//cout <<"Current iter " << iter << " Score after beta-theta " << newScore << endl;
					for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
					{
						cIter->second->clear();
						delete cIter->second;
					}
					subiter++;
					affectedVariables.clear();
					showid++;
					attemptedMoves++;
				}
				k++;	
				if((currGlobalScore-scorePremodule)<=convThreshold)
				{
					notConverged=false;
				}
				else
				{
					redefineModules_Global();
				}
				iter++;
				scorePremodule=currGlobalScore;
			}
			moduleiter++;
		}
		cout <<"Final Score " << currGlobalScore << endl;
		finalScores[f]=currGlobalScore;
	}
	return 0;
}

int
GraphLearner::getTopRegs(map<int,int>& topRegs)
{	
	map<int,int> regIDCnt;
	int max=0;
	for(int i=0;i<moveSet.size();i++)
	{
		MetaMove* m= moveSet[i];
		int r=m->getSrcVertex();
		if(regIDCnt.find(r)==regIDCnt.end())
		{
			regIDCnt[r]=1;
		}
		else
		{
			regIDCnt[r]=regIDCnt[r]+1;
		}
		if(regIDCnt[r]>max)
		{
			max=regIDCnt[r];
		}
	}
	for(map<int,int>::iterator rIter=regIDCnt.begin();rIter!=regIDCnt.end();rIter++)
	{
		if(rIter->second==max)
		{
			topRegs[rIter->first]=max;
		}
	}
	//cout << "Found " << topRegs.size() << " top regulators controlling " << max << " genes " << endl;
	regIDCnt.clear();
	return 0;
}

int
GraphLearner::loadLastUpdate()
{
	string& dirname=outLocMap[evMgrSet.begin()->first];
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/last.txt",dirname.c_str(),currFold);
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=varManager->getVariableSet();
	while(inFile.good())
	{	
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Variable* v=NULL;
		int titer=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				int vid=varManager->getVarID(tok);
				v=varSet[vid];
			}
			else
			{
				titer=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		variableStatus[v->getName()]=titer;
	}
	return 0;
}

int
GraphLearner::writeLastUpdate()
{
	string& dirname=outLocMap[evMgrSet.begin()->first];
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/last.txt",dirname.c_str(),currFold);
	ofstream outFile(aFName);
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* var=varSet[vIter->first];
		if(geneModuleID.find(var->getName())==geneModuleID.end())
		{
			continue;
		}
		outFile << var->getName() << "\t" << variableStatus[var->getName()] << endl;
	}
	outFile.close();
	return 0;
}

double
GraphLearner::loadInitPLLScore()
{
	double initScore=0;

	string& dirname=outLocMap[evMgrSet.begin()->first];
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/pll.txt",dirname.c_str(),currFold);
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=varManager->getVariableSet();
	INTDBLMAP* plls=new INTDBLMAP;
	currPLLMap[evMgrSet.begin()->first]=plls;
	while(inFile.good())
	{	
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Variable* v=NULL;
		double val=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				int vid=varManager->getVarID(tok);
				v=varSet[vid];
			}
			else
			{
				val=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		(*plls)[v->getID()] = val;
		initScore += val;
	}
	return initScore;
}

int
GraphLearner::writePLLScore()
{
	string& dirname=outLocMap[evMgrSet.begin()->first];
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/pll.txt",dirname.c_str(),currFold);
	ofstream outFile(aFName);
	VSET& varSet=varManager->getVariableSet();
	INTDBLMAP* plls = currPLLMap[evMgrSet.begin()->first];
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* var=varSet[vIter->first];
		if(geneModuleID.find(var->getName())==geneModuleID.end())
		{
			continue;
		}
		double newPLL_s=(*plls)[vIter->first];
		outFile << var->getName() << "\t" << newPLL_s << endl;
	}
	outFile.close();
	return 0;
}

double
GraphLearner::getInitPLLScore()
{
	double initScore=0;
	VSET& varSet=varManager->getVariableSet();
	//Initially we just sum up the marginal likelihoods from each condition
	INTDBLMAP* plls=new INTDBLMAP;
	currPLLMap[evMgrSet.begin()->first]=plls;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end())
		{
			continue;
		}
		Variable* var=varSet[vIter->first];
		if(geneModuleID.find(var->getName())==geneModuleID.end())
		{
			continue;
		}
		double newPLL_s=getNewPLLScore_Condition(-1,vIter->first,NULL);
		double priorScore=varNeighborhoodPrior[vIter->first];
		(*plls)[vIter->first]=newPLL_s+priorScore;
		initScore=initScore+(*plls)[vIter->first];
	}
	return initScore;
}



double
GraphLearner::getPLLScore()
{
	VSET& varSet=varManager->getVariableSet();
	double gScore=0;
	for(map<int,INTDBLMAP*>::iterator eIter=currPLLMap.begin();eIter!=currPLLMap.end();eIter++)
	{
		//EvidenceManager* evMgr=eIter->second;
		INTDBLMAP* plls=eIter->second;
		for(INTDBLMAP_ITER dIter=plls->begin();dIter!=plls->end();dIter++)
		{
			if(isnan(gScore) || isinf(gScore))
			{
				//cout << "Found nan/inf for variable " << dIter->first << endl;
			}
			double cscore = dIter->second;
			if(isnan(cscore) || isinf(cscore))
			{
				cout << "$$$ Found nan/inf for variable " << dIter->first << "(" << varSet[dIter->first]->getName() << ")" << endl;
			}
			else
			{
				gScore=gScore+cscore;
			}
		}
	}
	return gScore;
}


double 
GraphLearner::getPriorChange()
{
	double priorContrib=0;
	for(map<string,int>::iterator aIter=edgeUpdates.begin();aIter!=edgeUpdates.end();aIter++)
	{
		double edgeContrib=edgePresenceProb[aIter->first];
		priorContrib=priorContrib+(log(edgeContrib)-log(1-edgeContrib));
	}
	for(map<string,int>::iterator aIter=edgeRemoveUpdates.begin();aIter!=edgeRemoveUpdates.end();aIter++)
	{
		double edgeContrib=edgePresenceProb[aIter->first];
		priorContrib=priorContrib+(log(1-edgeContrib)-log(edgeContrib));
	}
	return priorContrib;
}

double 
GraphLearner::getInitPrior()
{
	double graphPrior=0;
	double edgePresence=1/(1+exp(-1*beta1));
	for(map<string,double>::iterator aIter=edgePresenceProb.begin();aIter!=edgePresenceProb.end();aIter++)
	{
		//graphPrior=graphPrior+log(1-edgePresence);
		graphPrior=graphPrior+log(1-aIter->second);
		if(isinf(graphPrior)|| isnan(graphPrior))
		{
			//cout <<"Graph prior is "<< graphPrior << " after " << aIter->first << " for " << aIter->second << endl;
		}
	}
	return graphPrior;
}

double
GraphLearner::getPLLScore_Datapoint(int datapoint)
{
	double pll=0;
	double wt=0;
	VSET& varSet=varManager->getVariableSet();
	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		EvidenceManager* evMgr=eIter->second;
		EMAP* evidMap=evMgr->getEvidenceAt(datapoint);
		for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
		{
			int vId=vIter->first;
			double cll=0;
			for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
			{
				INTINTMAP* cset=csIter->second;
				if((*cset)[eIter->first]==0)
				{
					continue;
				}
				int tempid=csIter->first;
				FactorGraph* fg=fgGraphSet[csIter->first];
				SlimFactor* sFactor=fg->getFactorAt(vId);
				Potential* sPot=sFactor->potFunc;
				double pval=sPot->getCondPotValueFor(evidMap);
				if(pval<1e-50)
				{
					pval=1e-50;
				}
				cll=cll+(pval*wt);
			}
			pll=pll+log(cll);
		}
	}

	return pll;
}

int
GraphLearner::clearFoldSpecData()
{
	//Clear existing graphs
	for(map<int,FactorGraph*>::iterator fIter=fgGraphSet.begin();fIter!=fgGraphSet.end();fIter++)
	{
		delete fIter->second;
	}
	fgGraphSet.clear();
	for(map<string,INTINTMAP*>::iterator eIter=edgeConditionMap.begin();eIter!=edgeConditionMap.end();eIter++)
	{
		eIter->second->clear();
		delete eIter->second;
	}
	edgeConditionMap.clear();

	for(map<int,INTDBLMAP*>::iterator plIter=currPLLMap.begin();plIter!=currPLLMap.end();plIter++)
	{
		plIter->second->clear();
		delete plIter->second;
	}
	currPLLMap.clear();
	edgeUpdates.clear();
	edgeRemoveUpdates.clear();
	return 0;
}

int
GraphLearner::initEdgeSet(bool validation)
{
	if(condsetMap.size()==0)
	{
		initCondsetMap_Nopool();
	}
	//Create a factorgraph for conditionset combination
	//It does not matter which factor manager we use now
	//for(map<int,INTINTMAP*>::iterator cIter=condsetMap.begin();cIter!=condsetMap.end();cIter++)
	for(map<int,FactorManager*>::iterator fIter=fgMgrSet.begin();fIter!=fgMgrSet.end();fIter++)
	{
		FactorManager* fMgr=fIter->second;
		FactorGraph* condspecGraph=fMgr->createInitialFactorGraph();
		fgGraphSet[fIter->first]=condspecGraph;
	}

	map<string,int> testedEdges;
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++)
	{
		Variable* u=varSet[uIter->first];
		if((restrictedVarList.size()>0) && (restrictedVarList.find(u->getName())==restrictedVarList.end()))
		{
			continue;
		}

		for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
		{
			if(uIter->first==vIter->first)
			{
				continue;
			}
			Variable* v=varSet[vIter->first];
			if(geneModuleID.find(v->getName())==geneModuleID.end())
			{	
				continue;
			}
			string edgeKey;
			//This is going to be a directed graph
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());
			if(strcmp(edgeKey.c_str(),"FBgn0013263\tFBgn0004170")==0)
			{
				//cout <<"Stop here" << endl;
			}
			if(strcmp(edgeKey.c_str(),"YIL118W\tYBR078W")==0)
			{
				//cout << "*************" << "I am adding this edge!\n" << endl;
			}
			INTINTMAP* conditionSet=new INTINTMAP;
			edgeConditionMap[edgeKey]=conditionSet;
			for(map<int,INTINTMAP*>::iterator cIter=condsetMap.begin();cIter!=condsetMap.end();cIter++)
			{
				(*conditionSet)[cIter->first]=0;
			}
			double initPrior=getEdgePrior(uIter->first,vIter->first);
			initPrior=1/(1+exp(-1*initPrior));
			if(initPrior<1e-6)
			{
				initPrior=1e-6;
			}
			if(initPrior==1)
			{
				initPrior=1-1e-6;
			}
			edgePresenceProb[edgeKey]=initPrior;
			if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end())
			{
				varNeighborhoodPrior[vIter->first]=log(1-initPrior);
			}
			else
			{
				varNeighborhoodPrior[vIter->first]=varNeighborhoodPrior[vIter->first]+log(1-initPrior);
			}
		}
	}
	cout <<"Restricted varlist size: " << restrictedVarList.size() << endl;
	int n=varSet.size();
	int r=restrictedVarList.size();
	int expEdgeCnt=((r*(r-1))/2) + (r*(n-r)) ;
	cout <<"Inited " << edgeConditionMap.size() << " edges. Expected " << expEdgeCnt << endl;
	testedEdges.clear();	
	for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
	{
		FactorGraph* condspecGraph=gIter->second;
		PotentialManager* potMgr=potMgrSet[gIter->first];
		//Init the potentials
		for(int f=0;f<condspecGraph->getFactorCnt();f++)
		{
			SlimFactor* sFactor=condspecGraph->getFactorAt(f);
			if (sFactor->potFunc != NULL)
			{
				delete sFactor->potFunc;
			}
			sFactor->potFunc=new Potential;
			sFactor->potFunc->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
			sFactor->potFunc->potZeroInit();
			potMgr->populatePotential(sFactor->potFunc,random);
			sFactor->potFunc->initMBCovMean();
		}
	}
	return 0;
}

int
GraphLearner::initCondsetMap_Nopool()
{
	int ind=0;
	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		INTINTMAP* cset=new INTINTMAP;

		for(map<int,EvidenceManager*>::iterator dIter=evMgrSet.begin();dIter!=evMgrSet.end();dIter++)
		{
			if(eIter==dIter)
			{
				(*cset)[dIter->first]=1;
			}
			else
			{
				(*cset)[dIter->first]=0;
			}
		}
		int currind=(int)pow(2.0,ind);
		condsetMap[currind]=cset;
		string condKey;
		genCondSetKey(*cset,condKey);
		condsetKeyIDMap[condKey]=currind;
		condsetIDKeyMap[currind]=condKey;
		ind=ind+1;
	}
	return 0;
}

INTINTMAP*
GraphLearner::getConditionSet(int cind)
{
	if(condsetMap.find(cind)==condsetMap.end())
	{
		cout <<"Did not find any condition sets associated with " << cind << endl;
		exit(0);
	}
	return condsetMap[cind];
}

int
GraphLearner::getPredictionError_Holdout(int foldid)
{
	VSET& varSet=varManager->getVariableSet();
	char foldoutDirName[1024];
	char aFName[1024];
	string& dirname=outLocMap[evMgrSet.begin()->first];
	sprintf(foldoutDirName,"%s/fold%d",dirname.c_str(),foldid);
	sprintf(aFName,"%s/prediction.txt",foldoutDirName);
	ofstream pFile(aFName);
	int totalEvid=holdoutEvMgr->getNumberOfEvidences();
	map<int,double> varPLL;
	for(int d=0;d<totalEvid;d++)
	{
		//for each gc, get the expected value of this datapoint
		EMAP* evidMap=holdoutEvMgr->getEvidenceAt(d);
		for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
		{
			FactorGraph* fg=fgGraphSet[csIter->first];
			for(map<string,int>::iterator vIter=reqdTargetList.begin();vIter!=reqdTargetList.end();vIter++)
			{
				if(strcmp(vIter->first.c_str(),"FBgn0004465")==0)
				{
					//cout <<"Stop here"<< endl;
				}
				int vId=varManager->getVarID(vIter->first.c_str());
				if(vId==-1)
				{
					continue;
				}
				Variable* v=varSet[vId];
				double cll=0;
				SlimFactor* sFactor=fg->getFactorAt(vId);
				Potential* sPot=sFactor->potFunc;
				if(sPot==NULL)
				{
					//cout <<"Found null for factor="<< sFactor->fId
					//	<< "variable=" <<varSet[sFactor->fId]->getName() << endl;
				}
				if(evidMap->find(vId)==evidMap->end())
				{
					//cout <<"Skipping " << vIter->first << endl;
					continue;
				}
				double pval=sPot->getCondPotValueFor(evidMap);
				if(pval<1e-50)
				{
					pval=1e-50;
				}
				if(isinf(pval) || isnan(pval))
				{
					//cout <<"Stop here. Found nan/inf for " << vIter->first << " dtpt "<< d << " cset " << csIter->first << endl;
				}
				cll=log(pval);
				if(varPLL.find(vId)==varPLL.end())
				{
					varPLL[vId]=cll;
				}
				else
				{
					varPLL[vId]=varPLL[vId]+cll;
				}
			}
		}
	}
	pFile <<"GeneName";
	for(int i=0;i<totalEvid;i++)
	{
		pFile <<"\t" <<i;
	}
	for(int i=0;i<totalEvid;i++)
	{
		pFile <<"\t" <<i;
	}
	pFile <<"\t0";
	for (int i=0;i<totalEvid;i++)
	{
		pFile <<"\t" << i;
	}
	pFile << endl;
	for(map<int,double>::iterator vIter=varPLL.begin();vIter!=varPLL.end();vIter++)
	{
		int vId=vIter->first;
		Variable* var=varSet[vId];
		
		pFile <<var->getName();
		int regulatorCnt=0;
		map<int,int> regsExpressed;
		//First the predicted time course
		for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
		{
			FactorGraph* fg=fgGraphSet[csIter->first];
			SlimFactor* sFactor=fg->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			regulatorCnt=sFactor->mergedMB.size();
			for(int i=0;i<totalEvid;i++)
			{
				EMAP* evidMap=holdoutEvMgr->getEvidenceAt(i);
				int predfrom=0;
				double predval=sPot->predictSample(evidMap,predfrom);
				pFile <<"\t" << predval;
				regsExpressed[i]=predfrom;
			}
		}
		//Then the true time course
		for(int i=0;i<totalEvid;i++)
		{
			EMAP* evidMap=holdoutEvMgr->getEvidenceAt(i);
			Evidence* evid=(*evidMap)[vId];
			pFile <<"\t" <<evid->getEvidVal();
		}
		pFile <<"\t" << regulatorCnt;
		for(int i=0;i<totalEvid;i++)
		{
			pFile <<"\t" << regsExpressed[i];
		}
		pFile <<endl;
	}
	pFile.close();
	varPLL.clear();
	return 0;
}


int
GraphLearner::getPredictionError_CrossValid(int foldid)
{
	VSET& varSet=varManager->getVariableSet();
	char foldoutDirName[1024];
	char aFName[1024];
	string& dirname=outLocMap[evMgrSet.begin()->first];
	sprintf(foldoutDirName,"%s/fold%d",dirname.c_str(),foldid);
	/*
	sprintf(aFName,"%s/predictionerr_cv.txt",foldoutDirName);
	ofstream oFile(aFName);
	sprintf(aFName,"%s/prediction_cv.txt",foldoutDirName);
	ofstream pFile(aFName);
	*/
	EvidenceManager* evMgr=evMgrSet.begin()->second;
	INTINTMAP& testSet=evMgr->getTestSet();
	map<int,double> varPLL;
	for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
	{
		//for each gc, get the expected value of this datapoint
		EMAP* evidMap=evMgr->getEvidenceAt(dIter->first);
		for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
		{
			FactorGraph* fg=fgGraphSet[csIter->first];
			//for(map<string,int>::iterator vIter=reqdTargetList.begin();vIter!=reqdTargetList.end();vIter++)
			for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
			{
				int vId=varManager->getVarID(vIter->first.c_str());
				if(vId==-1)
				{
					continue;
				}
				Variable* v=varSet[vId];
				double cll=0;
				SlimFactor* sFactor=fg->getFactorAt(vId);
				Potential* sPot=sFactor->potFunc;
				if(sPot==NULL)
				{
					//cout <<"Found null for factor="<< sFactor->fId
					//	<< " variable=" <<varSet[sFactor->fId]->getName() << endl;
				}
				double pval=sPot->getCondPotValueFor(evidMap);
				if(pval<1e-50)
				{
					pval=1e-50;
				}
				if(isinf(pval) || isnan(pval))
				{
					//cout <<"Stop here. Found nan/inf for " << vIter->first << " dtpt "<< dIter->first << " cset " << csIter->first << endl;
				}
				cll=log(pval);
				if(varPLL.find(vId)==varPLL.end())
				{
					varPLL[vId]=cll;
				}
				else
				{
					varPLL[vId]=varPLL[vId]+cll;
				}
			}
		}
	}
	/*
	for(map<int,double>::iterator pIter=varPLL.begin();pIter!=varPLL.end();pIter++)
	{
		oFile << varSet[pIter->first]->getName() << "\t" << pIter->second << endl;
	}
	pFile << "\tRMSE\tNormRMSE\tCoeff_Det_aka_R^2\tCC"<< endl;
	*/
	Distance d;
	vector<double> truevect;
	vector<double> predvect;
	for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
	{
		int vId=varManager->getVarID(vIter->first.c_str());
		if(vId==-1)
		{
			continue;
		}
		//pFile <<vIter->first;
		double error=0;
		double norm=0;
		double maxval=-100000;
		double minval=1000000;
		double totalvar=0;
		double truemean=0;
		truevect.clear();
		predvect.clear();
		for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
		{
			FactorGraph* fg=fgGraphSet[csIter->first];
			SlimFactor* sFactor=fg->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
			{
				EMAP* evidMap=evMgr->getEvidenceAt(dIter->first);
				Evidence* evid=(*evidMap)[vId];
				double trueval=evid->getEvidVal();
				truemean=truemean+trueval;
				truevect.push_back(trueval);
			}
		}
		truemean=truemean/((double)testSet.size());
		//First the predicted time course
		for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
		{
			FactorGraph* fg=fgGraphSet[csIter->first];
			SlimFactor* sFactor=fg->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
			{
				EMAP* evidMap=evMgr->getEvidenceAt(dIter->first);
				double predval=sPot->predictSample(evidMap);
				Evidence* evid=(*evidMap)[vId];
				double trueval=evid->getEvidVal();
				totalvar=totalvar+((trueval-truemean)*(trueval-truemean));
				//also called residuals
				error=error+((predval-trueval)*(predval-trueval));
				predvect.push_back(predval);
				//norm=norm+(trueval*trueval);
				norm=norm+1;
				if(trueval>maxval)
				{
					maxval=trueval;
				}
				if(trueval<minval)
				{
					minval=trueval;
				}
			}
		}
		//Then the true time course
		/*for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(dIter->first);
			Evidence* evid=(*evidMap)[vId];
			//pFile <<"\t" <<evid->getEvidVal();
		}*/
		double coeff_det=1-(error/totalvar);
		error=error/norm;
		error=sqrt(error);
		double nmsd=error/(maxval-minval);
		double cc=d.computeCC(truevect,predvect);
		//pFile<< "\t"<<error <<"\t" << nmsd << "\t" << coeff_det <<"\t" << cc <<endl;
	}
	//oFile.close();
	//pFile.close();
	varPLL.clear();
	return 0;
}

int
GraphLearner::collectMoves(int currK)
{
	VSET& varSet=varManager->getVariableSet();
	for(int i=0;i<moveSet.size();i++)
	{
		delete moveSet[i];
	}
	moveSet.clear();
	map<string,int> testedEdges;
	double step=1.0/(double)geneModuleID.size();
	double rVal=gsl_ran_flat(rnd,0,1);
	int rind=(int)(rVal/step);
	if(rind==varSet.size())
	{
		rind=varSet.size()-1;
	}
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* v=vIter->second;
		if(geneModuleID.find(v->getName())==geneModuleID.end())
		{
			continue;
		}
		
		int moduleID=geneModuleID[v->getName()];
		map<string,int>* moduleMembers=moduleGeneSet[moduleID];
		double bestTargetScore=0;
		double bestScoreImprovement=0;
		int bestCondSetInd=-1;
		Variable* bestu=NULL;
		map<string,INTDBLMAP*> impEdgeScoreImprovement;
		map<int,double> bestConditionImprovement;
		map<int,double> bestConditionScore;
		map<int,double> bestCondSpecPLL;
		Potential* bestPot=NULL;
		bool bestUpdate=false;
		for(map<string,int>::iterator uIter=restrictedVarList.begin();uIter!=restrictedVarList.end();uIter++)
		{
			int regID=varManager->getVarID(uIter->first.c_str());
			if(regID==-1)
			{
				continue;
			}
			if(vIter->first==regID)
			{
				continue;
			}
			Variable* u=varSet[regID];
			string edgeKey;
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());
			if(strcmp(edgeKey.c_str(),"YNL236W\tYPL188W")==0 || strcmp(edgeKey.c_str(),"YER068W\tYPR161C")==0)
			{
				//cout <<"Stop here" << endl;
			}
			if(testedEdges.find(edgeKey)!=testedEdges.end())
			{
				continue;
			}
			testedEdges[edgeKey]=0;
			//Generate next condition assignments
			if(edgeConditionMap.find(edgeKey)==edgeConditionMap.end())
			{
				cout <<"No edge " << edgeKey.c_str() << " u " << u->getID() << " v " << v->getID()<< endl;
				exit(0);
			}
			INTINTMAP* conditionSet=edgeConditionMap[edgeKey];
			map<int,double> conditionImprovement;
			map<int,double> conditionScore;
			for(INTINTMAP_ITER cIter=conditionSet->begin();cIter!=conditionSet->end();cIter++)
			{
				double scoreImprovement=0;
				double newTargetScore=0;
				Potential* aPot=NULL;
				INTDBLMAP pll;
				INTINTMAP* cset=getConditionSet(cIter->first);
				if(cIter->second==1)
				{
					pll.clear();
					continue;
				}
				if(!checkMBSize(cIter->first,regID,vIter->first,currK))
				{
					pll.clear();
					continue;
				}
				getNewPLLScore(cIter->first,*cset,u,v,edgeKey,newTargetScore,scoreImprovement,&aPot);
				double moduleWideScoreImprovement=getModuleWideScoreImprovement(cIter->first,*cset,u,v,moduleMembers);
				conditionImprovement[cIter->first]=scoreImprovement;
				conditionScore[cIter->first]=newTargetScore;
				if(scoreImprovement<0)
				{
					pll.clear();
					continue;
				}
				if((bestCondSetInd==-1) || (bestScoreImprovement<(scoreImprovement+moduleWideScoreImprovement)))
				{
					bestu=u;
					bestTargetScore=newTargetScore;
					bestScoreImprovement=scoreImprovement;
					bestCondSetInd=cIter->first;
					bestUpdate=true;
					bestCondSpecPLL.clear();
					for(INTDBLMAP_ITER dIter=pll.begin();dIter!=pll.end();dIter++)
					{
						bestCondSpecPLL[dIter->first]=dIter->second;
					}
					if(bestPot!=NULL)
					{
						delete bestPot;
					}
					bestPot=aPot;
				}
				else
				{
					delete aPot;
				}
				pll.clear();
			}
			if(bestUpdate)
			{
				bestConditionImprovement.clear();
				for(INTDBLMAP_ITER cvIter=conditionImprovement.begin();cvIter!=conditionImprovement.end();cvIter++)
				{
					bestConditionImprovement[cvIter->first]=cvIter->second;
					bestConditionScore[cvIter->first]=conditionScore[cvIter->first];
				}
			}
			conditionImprovement.clear();
			conditionScore.clear();
			bestUpdate=false;
		}
		//At this stage we have the best condition set for the pair {u,bestv}.
		if((bestu==NULL) || (bestScoreImprovement<=0))
		{
			continue;
		}
		MetaMove* move=new MetaMove;
		move->setSrcVertex(bestu->getID());
		move->setConditionSetInd(bestCondSetInd);
		for(INTDBLMAP_ITER dIter=bestCondSpecPLL.begin();dIter!=bestCondSpecPLL.end();dIter++)
		{
			move->pll[dIter->first]=dIter->second;
		}
		//cout <<"CurrK: "<< currK << " Found edge for " << bestCondSetInd <<" "  << bestu->getName().c_str() << "<->" << v->getName() << " score deltas";
		for(INTDBLMAP_ITER cvIter=bestConditionImprovement.begin();cvIter!=bestConditionImprovement.end();cvIter++)
		{
			//cout << " "<<cvIter->first<<"=" << cvIter->second;
		}
		//cout << " scores" ;
		for(INTDBLMAP_ITER cvIter=bestConditionScore.begin();cvIter!=bestConditionScore.end();cvIter++)
		{
			//cout << " "<<cvIter->first<<"=" << cvIter->second;
		}
		/*for(map<string,double>::iterator kIter=knownRegulatorScore.begin();kIter!=knownRegulatorScore.end();kIter++)
		{
			cout <<" " << kIter->first<<"="<< kIter->second;
		}*/
		//cout << endl;
		bestCondSpecPLL.clear();
		bestConditionScore.clear();
		move->setTargetVertex(v->getID());
		move->setTargetMBScore(bestTargetScore);
		move->setScoreImprovement(bestScoreImprovement);
		move->setDestPot(bestPot);
		moveSet.push_back(move);
	}
	testedEdges.clear();
	return 0;
}

int
GraphLearner::collectRemoveMoves(int currK,int rind)
{
	VSET& varSet=varManager->getVariableSet();
	for(int i=0;i<removeMoveSet.size();i++)
	{
		delete removeMoveSet[i];
	}
	removeMoveSet.clear();
	map<string,int> testedEdges;
	int vID=idVidMap[rind];
	VSET_ITER vIter=varSet.find(vID);
	if(vIter==varSet.end())
	{
		return 0;
	}
	//for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	//{
		Variable* v=vIter->second;
		if(geneModuleID.find(v->getName())==geneModuleID.end())
		{
			return 0;
			//continue;
		}
		int moduleID=geneModuleID[v->getName()];
		map<string,int>* moduleMembers=moduleGeneSet[moduleID];
		double bestTargetScore=0;
		double bestScoreImprovement=0;
		int bestCondSetInd=-1;
		Variable* bestu=NULL;
		map<int,double> bestCondSpecPLL;
		Potential* bestPot=NULL;
		bool bestUpdate=false;
		INTINTMAP_ITER myCIter=condsetMap.begin()->second->begin();
		PotentialManager* potMgr=potMgrSet[myCIter->first];
		FactorGraph* fg=fgGraphSet[myCIter->first];
		SlimFactor* dFactor=fg->getFactorAt(v->getID());
		vector<int> regvec;
		for(INTINTMAP_ITER uIter=dFactor->mergedMB.begin();uIter!=dFactor->mergedMB.end();uIter++)
		{
			int regID=uIter->first;
			regvec.push_back(regID);
		}
		if (regvec.size() == 1)
		{
			return 0;
		}
		for(int uIter=0;uIter<regvec.size();uIter++)
		{
			int regID=regvec[uIter];
			if(vIter->first==regID)
			{
				continue;
			}
			Variable* u=varSet[regID];
			string edgeKey;
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());
			if (removedEdges.find(edgeKey) != removedEdges.end())
			{
				continue;
			}
			if (addedEdges.find(edgeKey) != addedEdges.end())
			{
				continue;
			}
			if(testedEdges.find(edgeKey)!=testedEdges.end())
			{
				continue;
			}
			testedEdges[edgeKey]=0;
			//Generate next condition assignments
			if(edgeConditionMap.find(edgeKey)==edgeConditionMap.end())
			{
				cout <<"(collect remove) No edge " << edgeKey.c_str() << " u (" << u->getName() << ") " << u->getID() ;
				cout << " v (" << v->getName() << ") " << v->getID()<< endl;
				exit(0);
			}
			INTINTMAP* conditionSet=edgeConditionMap[edgeKey];
			for(INTINTMAP_ITER cIter=conditionSet->begin();cIter!=conditionSet->end();cIter++)
			{
				double scoreImprovement=0;
				double newTargetScore=0;
				Potential* aPot=NULL;
				INTDBLMAP pll;
				INTINTMAP* cset=getConditionSet(cIter->first);
				if(cIter->second==1)
				{
					pll.clear();
					continue;
				}
				getPLLScoreRemoveOne(cIter->first,*cset,u,v,edgeKey,newTargetScore,scoreImprovement,&aPot);
				//double moduleWideScoreImprovement=getModuleWideScoreImprovement(cIter->first,*cset,u,v,moduleMembers);
				double moduleWideScoreImprovement=0;
				if(scoreImprovement<0)
				{
					pll.clear();
					continue;
				}
				if((bestCondSetInd==-1) || (bestScoreImprovement<(scoreImprovement+moduleWideScoreImprovement)))
				{
					bestu=u;
					bestTargetScore=newTargetScore;
					bestScoreImprovement=scoreImprovement;
					bestCondSetInd=cIter->first;
					bestUpdate=true;
					bestCondSpecPLL.clear();
					for(INTDBLMAP_ITER dIter=pll.begin();dIter!=pll.end();dIter++)
					{
						bestCondSpecPLL[dIter->first]=dIter->second;
					}
					if(bestPot!=NULL)
					{
						delete bestPot;
					}
					bestPot=aPot;
				}
				else
				{
					delete aPot;
				}
				pll.clear();
			}
			bestUpdate=false;
		}
		//At this stage we have the best condition set for the pair {u,bestv}.
		if((bestu==NULL) || (bestScoreImprovement<=0))
		{
			testedEdges.clear();
			return 0;
		//	continue;
		}
		MetaMove* move=new MetaMove;
		move->setSrcVertex(bestu->getID());
		move->setConditionSetInd(bestCondSetInd);
		for(INTDBLMAP_ITER dIter=bestCondSpecPLL.begin();dIter!=bestCondSpecPLL.end();dIter++)
		{
			move->pll[dIter->first]=dIter->second;
		}
		bestCondSpecPLL.clear();
		move->setTargetVertex(v->getID());
		move->setTargetMBScore(bestTargetScore);
		move->setScoreImprovement(bestScoreImprovement);
		move->setDestPot(bestPot);
		removeMoveSet.push_back(move);
	//}
	testedEdges.clear();
	return 0;
}

int
GraphLearner::collectMoves(int currK,int rind)
{
	VSET& varSet=varManager->getVariableSet();
	for(int i=0;i<moveSet.size();i++)
	{
		delete moveSet[i];
	}
	moveSet.clear();
	map<string,int> testedEdges;
	int vID=idVidMap[rind];
	VSET_ITER vIter=varSet.find(vID);
	if(vIter==varSet.end())
	{
		return 0;
	}
	//for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	//{
		Variable* v=vIter->second;
		if(geneModuleID.find(v->getName())==geneModuleID.end())
		{
			return 0;
			//continue;
		}
		int moduleID=geneModuleID[v->getName()];
		map<string,int>* moduleMembers=moduleGeneSet[moduleID];
		double bestTargetScore=0;
		double bestScoreImprovement=0;
		int bestCondSetInd=-1;
		Variable* bestu=NULL;
		map<int,double> bestCondSpecPLL;
		Potential* bestPot=NULL;
		bool bestUpdate=false;
		for(map<string,int>::iterator uIter=restrictedVarList.begin();uIter!=restrictedVarList.end();uIter++)
		{
			int regID=varManager->getVarID(uIter->first.c_str());
			if(regID==-1)
			{
				continue;
			}
			
			if(vIter->first==regID)
			{
				continue;
			}
			Variable* u=varSet[regID];
			/*if(priorGraph_ChIP.find(u->getName())==priorGraph_ChIP.end())
			{
				continue;
			}
			map<string,double>* tgts=priorGraph_ChIP[u->getName()];
			if(tgts->find(v->getName())==tgts->end())
			{
				continue;
			}*/
			string edgeKey;
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());

			if (addedEdges.find(edgeKey) != addedEdges.end())
			{
				continue;
			}
			if (removedEdges.find(edgeKey) != removedEdges.end())
			{
				continue;
			}
			if(testedEdges.find(edgeKey)!=testedEdges.end())
			{
				continue;
			}
			if(resumeEdges.find(edgeKey)!=resumeEdges.end())
			{
				continue;
			}
			testedEdges[edgeKey]=0;
			//Generate next condition assignments
			if(edgeConditionMap.find(edgeKey)==edgeConditionMap.end())
			{
				cout <<"(collect move) No edge " << edgeKey.c_str() << " u " << u->getID() << " v " << v->getID()<< endl;
				exit(0);
			}
			INTINTMAP* conditionSet=edgeConditionMap[edgeKey];
			for(INTINTMAP_ITER cIter=conditionSet->begin();cIter!=conditionSet->end();cIter++)
			{
				double scoreImprovement=0;
				double newTargetScore=0;
				Potential* aPot=NULL;
				INTDBLMAP pll;
				INTINTMAP* cset=getConditionSet(cIter->first);
				if(cIter->second==1)
				{
					pll.clear();
					continue;
				}
				if(!checkMBSize(cIter->first,regID,vIter->first,currK))
				{
					pll.clear();
					continue;
				}
				getNewPLLScore(cIter->first,*cset,u,v,edgeKey,newTargetScore,scoreImprovement,&aPot);
				//double moduleWideScoreImprovement=getModuleWideScoreImprovement(cIter->first,*cset,u,v,moduleMembers);
				if(scoreImprovement<0)
				{
					pll.clear();
					continue;
				}
				//if((bestCondSetInd==-1) || (bestScoreImprovement<(scoreImprovement+moduleWideScoreImprovement)))
				if((bestCondSetInd==-1) || (bestScoreImprovement<(scoreImprovement)))
				{
					bestu=u;
					bestTargetScore=newTargetScore;
					bestScoreImprovement=scoreImprovement;
					bestCondSetInd=cIter->first;
					bestUpdate=true;
					bestCondSpecPLL.clear();
					for(INTDBLMAP_ITER dIter=pll.begin();dIter!=pll.end();dIter++)
					{
						bestCondSpecPLL[dIter->first]=dIter->second;
					}
					if(bestPot!=NULL)
					{
						delete bestPot;
					}
					bestPot=aPot;
					if(strcmp(v->getName().c_str(),"YHR108W")==0)
					{
						//cout <<"Stop here" << endl;
					}
				}
				else
				{
					delete aPot;
				}
				pll.clear();
			}
			bestUpdate=false;
		}
		//At this stage we have the best condition set for the pair {u,bestv}.
		if((bestu==NULL) || (bestScoreImprovement<=0))
		{
			testedEdges.clear();
			return 0;
		//	continue;
		}
		MetaMove* move=new MetaMove;
		move->setSrcVertex(bestu->getID());
		move->setConditionSetInd(bestCondSetInd);
		for(INTDBLMAP_ITER dIter=bestCondSpecPLL.begin();dIter!=bestCondSpecPLL.end();dIter++)
		{
			move->pll[dIter->first]=dIter->second;
		}
		/*cout <<"CurrK: "<< currK << " Found edge for " << bestCondSetInd <<" "  << bestu->getName().c_str() << "<->" << v->getName() << " score deltas";
		cout << " scores" ;
		for(map<string,double>::iterator kIter=knownRegulatorScore.begin();kIter!=knownRegulatorScore.end();kIter++)
		{
			cout <<" " << kIter->first<<"="<< kIter->second;
		}
		cout << endl; */
		bestCondSpecPLL.clear();
		move->setTargetVertex(v->getID());
		move->setTargetMBScore(bestTargetScore);
		move->setScoreImprovement(bestScoreImprovement);
		move->setDestPot(bestPot);
		moveSet.push_back(move);
	//}
	testedEdges.clear();
	return 0;
}
double
GraphLearner::getModuleWideScoreImprovement(int cid,INTINTMAP& conditionSet,Variable* u, Variable* v,map<string,int>* moduleMembers)
{
	return 0;
	VSET& varSet=varManager->getVariableSet();
	int vID=v->getID();
	int uID=u->getID();
	double moduleScore=0;
	double neighborCnt=0;
	FactorGraph* currGraph=fgGraphSet.begin()->second;
	for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
	{
		int nID=varManager->getVarID(mIter->first.c_str());
		if(nID==vID)
		{
			continue;
		}

		if(nID==uID)
		{
			continue;
		}
		SlimFactor* sFactor=currGraph->getFactorAt(nID);
		if(sFactor->mergedMB.find(uID)!=sFactor->mergedMB.end())
		{
			continue;
		}
		double scoreImprovement;
		double mbScore;
		Potential* dPot=NULL;
		string edgeKey(u->getName());
		edgeKey.append("\t");
		edgeKey.append(mIter->first.c_str());
		getNewPLLScore(cid,conditionSet,u,varSet[nID],edgeKey,mbScore,scoreImprovement,&dPot);
		if(dPot!=NULL)
		{
			delete dPot;
		}
		moduleScore=moduleScore+scoreImprovement;
		neighborCnt++;
	}
	double score=moduleScore/neighborCnt;
	return score;
}

int
GraphLearner::getPLLScoreRemoveOne(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot)
{
	VSET& varSet=varManager->getVariableSet();
	map<int,Potential*> dPots;
	scoreImprovement=0;
	bool delFromD=true;
	double currPrior=varNeighborhoodPrior[v->getID()];
	double plus=0;
	double minus=0;
	SlimFactor* tfactor;
	for(INTINTMAP_ITER cIter=conditionSet.begin();cIter!=conditionSet.end();cIter++)
	{
		PotentialManager* potMgr=potMgrSet[cIter->first];
		FactorGraph* fg=fgGraphSet[cIter->first];
		SlimFactor* dFactor=fg->getFactorAt(v->getID());
		tfactor = dFactor;
		//If it is still there, remove it
		INTINTMAP_ITER dIter=dFactor->mergedMB.find(u->getID());
		if(dIter!=dFactor->mergedMB.end())
		{
			dFactor->mergedMB.erase(dIter);
		}
		Potential *dPot=new Potential;
		dPot->setAssocVariable(varSet[dFactor->fId],Potential::FACTOR);
		for(INTINTMAP_ITER mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++)
		{
			Variable* aVar=varSet[mIter->first];
			dPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			double eprior=getEdgePrior(mIter->first,v->getID());
			double moduleContrib=getModuleContribLogistic((string&)v->getName(),(string&)aVar->getName());
			double edgeProb=1/(1+exp(-1*(eprior+moduleContrib)));
			//double edgeProb=1/(1+exp(-1*(eprior)));
			double edgeProbOld=1/(1+exp(-1*(eprior)));
			minus=minus+log(1-edgeProbOld);
			plus=plus+log(edgeProb);
		}
		dPots[cIter->first]=dPot;
		dPot->potZeroInit();
		dPot->setCondBias(dFactor->potFunc->getCondBias());
		dPot->setCondVariance(dFactor->potFunc->getCondVariance());
		//dPot->setCondWeight(dFactor->potFunc->getCondWeight());
	}
	currPrior=currPrior+plus-minus;
	string condKey;
	genCondSetKey(conditionSet,condKey);
	PotentialManager* potMgr=potMgrSet[cid];
	double newPLL_d=0;
	*newdPot=dPots[cid];
	potMgr->populatePotential(*newdPot,random);
	(*newdPot)->initMBCovMean();
	INTDBLMAP& regWts=(*newdPot)->getCondWeight();
	//cout << "&&& testing rm: " << tfactor->mergedMB.size() << "," << regWts.size() << endl;
	for(map<int,Potential*>::iterator pIter=dPots.begin();pIter!=dPots.end();pIter++)
	{
		Potential* dPot=pIter->second;
		if((dPot->getCondVariance()<0) || (isnan(dPot->getCondVariance())) || (isinf(dPot->getCondVariance())))
		{
			scoreImprovement=-1;
		}
	}
	if(scoreImprovement!=-1)
	{
		//newPLL_d=getNewPLLScore_Condition(cid,v->getID(),*newdPot);
		newPLL_d=getNewPLLScore_Condition(cid,v->getID(),u->getID(),*newdPot);
		if (isnan(newPLL_d))
		{
			cout << "$$$ there is a NAN (when removing, before prior): " << edgeKey << endl;
		}
		//newPLL_d=getNewPLLScore_Condition_Tracetrick(cid,v->getID(),u->getID(),*newdPot);
		newPLL_d=newPLL_d+currPrior;
		int keyid=evMgrSet.begin()->first;
		INTDBLMAP* plls=currPLLMap[keyid];
		double oldPLL_d=(*plls)[v->getID()];
		double dImpr=newPLL_d-oldPLL_d;
		if(edgePresenceProb.find(edgeKey)==edgePresenceProb.end())
		{
			cout <<"No edge prior for " << edgeKey.c_str() << endl;
			exit(0);
		}
	//	double priorImpr=log(edgeProb)-log(1-edgeProb);
		//dImpr=dImpr+priorImpr;
		
		if(dImpr<=0)
		{
			scoreImprovement=-1;
		}
		else
		{
			if(!delFromD)
			{	
				//cout <<"Trying to add the same regulator " << u->getName() <<" to " << v->getName() << endl;
			}

			scoreImprovement=dImpr;
			//mbScore=newPLL_d + log(edgeProb);
			mbScore=newPLL_d;
			if (isnan(newPLL_d))
			{
				cout << "$$$ there is a NAN (when removing): " << edgeKey << endl;
			}
			for(map<int,Potential*>::iterator sIter=dPots.begin();sIter!=dPots.end();sIter++)
			{
				if(sIter->first==cid)
				{
					continue;
				}
				delete sIter->second;
			}
			dPots.clear();
		}
	}
	//Put U back
	for(INTINTMAP_ITER cIter=conditionSet.begin();cIter!=conditionSet.end();cIter++)
	{
		FactorGraph* fg=fgGraphSet[cIter->first];
		SlimFactor* dFactor=fg->getFactorAt(v->getID());
		dFactor->mergedMB[u->getID()]=0;
	}

	if(scoreImprovement<0)
	{
		for(map<int,Potential*>::iterator pIter=dPots.begin();pIter!=dPots.end();pIter++)
		{
			delete pIter->second;
		}
		dPots.clear();
		*newdPot=NULL;
	}
	return 0;
}

int
GraphLearner::getNewPLLScore(int cid, INTINTMAP& conditionSet, Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot)
{
	VSET& varSet=varManager->getVariableSet();
	map<int,Potential*> dPots;
	map<int,bool> dPotDels;
	scoreImprovement=0;
	bool delFromD=true;
	double currPrior=varNeighborhoodPrior[v->getID()];
	double plus=0;
	double minus=0;
	for(INTINTMAP_ITER cIter=conditionSet.begin();cIter!=conditionSet.end();cIter++)
	{
		PotentialManager* potMgr=potMgrSet[cIter->first];
		FactorGraph* fg=fgGraphSet[cIter->first];
		SlimFactor* dFactor=fg->getFactorAt(v->getID());
		//Pretend as if we were already adding dFactor into sFactor's MB
		if(cIter->second==1)
		{
			if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
			{
				delFromD=false;
			}
			dPotDels[cIter->first]=delFromD;
			dFactor->mergedMB[u->getID()]=0;
		}
		Potential *dPot=new Potential;
		dPot->setAssocVariable(varSet[dFactor->fId],Potential::FACTOR);
		for(INTINTMAP_ITER mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++)
		{
			Variable* aVar=varSet[mIter->first];
			dPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			double eprior=getEdgePrior(mIter->first,v->getID());
			//double moduleContrib=getModuleContribLogistic((string&)v->getName(),(string&)u->getName());
			double moduleContrib=getModuleContribLogistic((string&)v->getName(),(string&)aVar->getName());
			double edgeProb=1/(1+exp(-1*(eprior+moduleContrib)));
			//double edgeProb=1/(1+exp(-1*(eprior)));
			double edgeProbOld=1/(1+exp(-1*(eprior)));
			minus=minus+log(1-edgeProbOld);
			plus=plus+log(edgeProb);
		}
		dPots[cIter->first]=dPot;
		dPot->potZeroInit();
		dPot->setCondBias(dFactor->potFunc->getCondBias());
		dPot->setCondVariance(dFactor->potFunc->getCondVariance());
		dPot->setCondWeight(dFactor->potFunc->getCondWeight());
	}
	currPrior=currPrior+plus-minus;
	string condKey;
	genCondSetKey(conditionSet,condKey);
	PotentialManager* potMgr=potMgrSet[cid];
	double newPLL_d=0;
	*newdPot=dPots[cid];
	potMgr->populatePotential(*newdPot,random);
	(*newdPot)->initMBCovMean();
	for(map<int,Potential*>::iterator pIter=dPots.begin();pIter!=dPots.end();pIter++)
	{
		Potential* dPot=pIter->second;
		if((dPot->getCondVariance()<0) || (isnan(dPot->getCondVariance())) || (isinf(dPot->getCondVariance()))
		)
		{
			scoreImprovement=-1;
		}
	}
	if(scoreImprovement!=-1)
	{
		//newPLL_d=getNewPLLScore_Condition(cid,v->getID(),*newdPot);
		//newPLL_d=getNewPLLScore_Condition(cid,v->getID(),u->getID(),*newdPot);
		newPLL_d=getNewPLLScore_Condition_Tracetrick(cid,v->getID(),u->getID(),*newdPot);
		if (isnan(newPLL_d))
		{
			cout << "$$$ there is a NAN (when adding, before prior): " << edgeKey << endl;
		}
		newPLL_d=newPLL_d+currPrior;
		int keyid=evMgrSet.begin()->first;
		INTDBLMAP* plls=currPLLMap[keyid];
		double oldPLL_d=(*plls)[v->getID()];
		double dImpr=newPLL_d-oldPLL_d;
		if(edgePresenceProb.find(edgeKey)==edgePresenceProb.end())
		{
			cout <<"No edge prior for " << edgeKey.c_str() << endl;
			exit(0);
		}
	//	double priorImpr=log(edgeProb)-log(1-edgeProb);
		//dImpr=dImpr+priorImpr;
		
		if(dImpr<=0)
		{
			scoreImprovement=-1;
		}
		else
		{
			if(!delFromD)
			{	
				//cout <<"Trying to add the same regulator " << u->getName() <<" to " << v->getName() << endl;
			}

			scoreImprovement=dImpr;
			//mbScore=newPLL_d + log(edgeProb);
			mbScore=newPLL_d;
			if (isnan(newPLL_d))
			{
				cout << "$$$ there is a NAN (when adding): " << edgeKey << endl;
			}
			for(map<int,Potential*>::iterator sIter=dPots.begin();sIter!=dPots.end();sIter++)
			{
				if(sIter->first==cid)
				{
					continue;
				}
				delete sIter->second;
			}
			dPots.clear();
		}
	}
	for(map<int,bool>::iterator bIter=dPotDels.begin();bIter!=dPotDels.end();bIter++)
	{
		if(bIter->second==false)
		{
			continue;
		}
		bool delFromD=bIter->first;
		FactorGraph* fg=fgGraphSet[bIter->first];
		SlimFactor* dFactor=fg->getFactorAt(v->getID());
		INTINTMAP_ITER dIter=dFactor->mergedMB.find(u->getID());
		dFactor->mergedMB.erase(dIter);
	}
	if(scoreImprovement<0)
	{
		for(map<int,Potential*>::iterator pIter=dPots.begin();pIter!=dPots.end();pIter++)
		{
			delete pIter->second;
		}
		dPots.clear();
		*newdPot=NULL;
	}
	return 0;
}



double
GraphLearner::getNewPLLScore_Condition(int csetId, int vId, Potential* newPot)
{
	double pll=0; 
	//Need to fix this to be set automatically
	double paramCnt=0;
	int datasize=0;
	//get parameter prior
	double paramPrior=0;
	for(map<int,EvidenceManager*>::iterator dIter=evMgrSet.begin();dIter!=evMgrSet.end();dIter++)
	{
		EvidenceManager* evMgr=evMgrSet[dIter->first];
		INTINTMAP* tSet=NULL;
		tSet=&evMgr->getTrainingSet();
		datasize=datasize+tSet->size();
		for(INTINTMAP_ITER eIter=tSet->begin();eIter!=tSet->end();eIter++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
			//Go over all condition sets that include cInd
			double cll=0;
			for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
			{
				INTINTMAP* cset=csIter->second;
				if((*cset)[dIter->first]==0)
				{
					continue;
				}
				FactorGraph* fg=fgGraphSet[csIter->first];
				SlimFactor* sFactor=fg->getFactorAt(vId);
				Potential* sPot=sFactor->potFunc;
				if((csIter->first==csetId) && (newPot!=NULL))
				{
					sPot=newPot;
				}
				double pval=sPot->getCondPotValueFor(evidMap);
				if(isnan(pval))
				{
					//cout <<"Pval is nan for condition " << csIter->first << " datapoint " << eIter->first << endl;
				}
				if(pval<1e-50)
				{
					pval=1e-50;
				}
				cll=cll+pval;
				if((eIter==tSet->begin()) && (dIter==evMgrSet.begin()))
				{
					double vCnt=(double)sPot->getAssocVariables().size();
					paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
					//paramCnt=vCnt-1;
				}
			}
			pll=pll+log(cll);
		}
	}	
	if(newPot==NULL)
	{
		FactorGraph* fg=fgGraphSet.begin()->second;
		SlimFactor* sFactor=fg->getFactorAt(vId);
		Potential* sPot=sFactor->potFunc;
		double testll=sPot->computeLL_Tracetrick(datasize);
		//cout <<"TestLL: " << testll << " pll " << pll << endl; 
	}
	//pll=pll-(0.5*paramCnt*log(datasize));
	pll=pll-(lambda*paramCnt*log(datasize));
	
	return pll;

}


double
GraphLearner::getNewPLLScore_Condition(int csetId, int vId, int uId, Potential* newPot)
{
	double pll=0; 
	//Need to fix this to be set automatically
	double paramCnt=0;
	int datasize=0;
	//get parameter prior
	double paramPrior=0;
	for(map<int,EvidenceManager*>::iterator dIter=evMgrSet.begin();dIter!=evMgrSet.end();dIter++)
	{
		EvidenceManager* evMgr=evMgrSet[dIter->first];
		INTINTMAP* tSet=NULL;
		tSet=&evMgr->getTrainingSet();
		datasize=datasize+tSet->size();
		for(INTINTMAP_ITER eIter=tSet->begin();eIter!=tSet->end();eIter++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
			//Go over all condition sets that include cInd
			double cll=0;
			for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
			{
				INTINTMAP* cset=csIter->second;
				if((*cset)[dIter->first]==0)
				{
					continue;
				}
				FactorGraph* fg=fgGraphSet[csIter->first];
				SlimFactor* sFactor=fg->getFactorAt(vId);
				Potential* sPot=sFactor->potFunc;
				if((csIter->first==csetId) && (newPot!=NULL))
				{
					sPot=newPot;
				}
				INTDBLMAP& cacheMeans=sFactor->cachePartialMean;
				double pval=0;
				if(cacheMeans.size()==0)
				{
					pval=sPot->getCondPotValueFor(evidMap,eIter->first);
				}
				else
				{
					//pval=sPot->getCondPotValueFor(evidMap,uId,cacheMeans,eIter->first);
					//double checkpval=sPot->getCondPotValueFor(evidMap);
					pval=sPot->getCondPotValueFor(evidMap);
				}
				if(isnan(pval))
				{
					//cout <<"Pval is nan for condition " << csIter->first << " datapoint " << eIter->first << endl;
				}
				if(pval<1e-50)
				{
					pval=1e-50;
				}
				cll=cll+pval;
				if((eIter==tSet->begin()) && (dIter==evMgrSet.begin()))
				{
					double vCnt=(double)sPot->getAssocVariables().size();
					paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
					//paramCnt=vCnt-1;
				}
			}
			pll=pll+log(cll);
		}
	}
	//pll=pll-(0.5*paramCnt*log(datasize));
	pll=pll-(lambda*paramCnt*log(datasize));
	return pll;
}


double
GraphLearner::getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uId, Potential* newPot)
{
	double pll=0; 
	//Need to fix this to be set automatically
	//get parameter prior
	double paramPrior=0;
	Potential* parentPot=new Potential;
	VSET& vars=newPot->getAssocVariables();
	for(VSET_ITER vIter=vars.begin();vIter!=vars.end();vIter++)
	{
		if(vIter->first==vId)
		{
			continue;
		}
		Variable* aVar=vars[vIter->first];
		parentPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	parentPot->potZeroInit();
	PotentialManager* potMgr=potMgrSet.begin()->second;
	potMgr->populatePotential(parentPot,false);
	//parentPot->initMBCovMean();
	EvidenceManager* evMgr=evMgrSet.begin()->second;
	INTINTMAP* tSet=&evMgr->getTrainingSet();
	int datasize=tSet->size();
	double jointll1=newPot->computeLL_Tracetrick(datasize);
	double jointll2=parentPot->computeLL_Tracetrick(datasize);
	double vCnt=(double)newPot->getAssocVariables().size();
	double paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
	pll=jointll1-jointll2;
	pll=pll-(lambda*paramCnt*log(datasize));
	delete parentPot;
	return pll;
}

double 
GraphLearner::getEdgePrior(int tfID, int targetID)
{
	INTDBLMAP* regPriors=NULL;
	double prior=beta1;
	double fwt = 0;
	for (map<string,map<int,INTDBLMAP*>*>::iterator pItr=edgepriormap.begin(); pItr!=edgepriormap.end(); pItr++)
	{
		double eweight=0;
		double gbeta = 0;
		map<int,INTDBLMAP*>* edgeprior = pItr->second;
		if(edgeprior->find(targetID)!=edgeprior->end())
		{
			regPriors=(*edgeprior)[targetID];
			if(regPriors->find(tfID)!=regPriors->end())
			{
				eweight=(*regPriors)[tfID];
				gbeta = betamap[pItr->first];
				fwt = fwt + gbeta*eweight;
			}
		}
	}
	prior=beta1+fwt;
	//if(prior<1e-6)
	//{
	//	prior=1e-6;
	//}
	//if(prior==1)
	//{
	//	prior=1-1e-6;
	//}
	return prior;
}

int
GraphLearner::genCondSetKey(INTINTMAP& condSet, string& aKey)
{	
	char keypair[256];
	for(INTINTMAP_ITER cIter=condSet.begin();cIter!=condSet.end();cIter++)
	{
		sprintf(keypair,"-%d=%d",cIter->first,cIter->second);
		aKey.append(keypair);
	}
	return 0;
}

int 
GraphLearner::sortMoves()
{
	for(int m=0;m<moveSet.size();m++)
	{
		for(int n=m+1;n<moveSet.size();n++)
		{
			MetaMove* m1=moveSet[m];
			MetaMove* m2=moveSet[n];
			if(m1->getScoreImprovement()<m2->getScoreImprovement())
			{
				moveSet[m]=m2;
				moveSet[n]=m1;
			}
		}
	}
	return 0;
}

int
GraphLearner::makeRemoveMoves()
{
	//map<int,INTINTMAP*> affectedVariables;
	for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
	{
		int cind=gIter->first;
		INTINTMAP* csVars=new INTINTMAP;
		affectedVariables[cind]=csVars;
	}
	int successMove=0;
	for(int m=0;m<removeMoveSet.size();m++)
	{
		MetaMove* move=removeMoveSet[m];
		int pool=0;
		if(attemptRemoveMove(move,affectedVariables)==0)
		{
			successMove++;
		}
		else
		{
			Potential* apot=move->getDestPot();
			delete apot;
		}
	}
	return 0;
}

int
GraphLearner::makeMoves()
{
	//map<int,INTINTMAP*> affectedVariables;
	for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
	{
		int cind=gIter->first;
		INTINTMAP* csVars=new INTINTMAP;
		affectedVariables[cind]=csVars;
	}
	int successMove=0;
	double netScoreDelta=0;
	for(int m=0;m<moveSet.size();m++)
	{
		MetaMove* move=moveSet[m];
		int pool=0;
		INTINTMAP* cSet=getConditionSet(move->getConditionSetInd());
		if(attemptMove(move,affectedVariables)==0)
		{
			successMove++;
			netScoreDelta=netScoreDelta+move->getScoreImprovement();
		}
		else
		{
			Potential* apot=move->getDestPot();
			delete apot;
		}
	}
	//cout <<"Total successful moves " << successMove << " out of total " << moveSet.size() << " with net score improvement " << netScoreDelta<< endl;
	return 0;
}


int
GraphLearner::attemptRemoveMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
	VSET& varSet=varManager->getVariableSet();
	string edgeKey;
	Variable* u=varSet[move->getSrcVertex()];
	Variable* v=varSet[move->getTargetVertex()];
	edgeKey.append(u->getName().c_str());
	edgeKey.append("\t");
	edgeKey.append(v->getName().c_str());
	string tgtName = v->getName();
	string tfName  = u->getName();
	
	if(edgeConditionMap.find(edgeKey)==edgeConditionMap.end())
	{
		//cout << "### didnt remove, edgeConditionMap!" << endl;
		//cout <<"Edge " << edgeKey << " not found " << endl;
		return -1;
	}
	INTINTMAP* conditionSet=edgeConditionMap[edgeKey];
	int csetind=move->getConditionSetInd();
	if(affectedVars.find(1)==affectedVars.end())
	{
		cout << "No csvars for condition " <<csetind<< endl;
		exit(0);
	}
	INTINTMAP* csVars=affectedVars[1];
	if((csVars->find(move->getSrcVertex())!=csVars->end()) || (csVars->find(move->getTargetVertex())!=csVars->end()))
	{
		//cout << "### didnt remove, affectedVars!" << endl;
		return -1;
	}
	removedEdges[edgeKey]=0;
	//(*csVars)[move->getSrcVertex()]=0;
	(*csVars)[move->getTargetVertex()]=0;
	INTINTMAP* condset=condsetMap[csetind];
	for(INTINTMAP_ITER cIter=condset->begin();cIter!=condset->end();cIter++)
	{
		if(cIter->second==0)
		{
			continue;
		}
		//cout << "### Attempted remove!" << endl;
		FactorGraph* csGraph=fgGraphSet[cIter->first];
		SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
		INTINTMAP_ITER dItr = dFactor->mergedMB.find(move->getSrcVertex());
		if(dItr!=dFactor->mergedMB.end())
		{
			//cout << "### Succeeded in remove!" << endl;
			//cout << "### size before rm: " << dFactor->mergedMB.size() << "," << edgeKey << endl;
			dFactor->mergedMB.erase(dItr);
			//cout << "### size after rm: " << dFactor->mergedMB.size() << "," << edgeKey << endl;
			//cout << "### pot size before rm: " << regWts.size() << "," << edgeKey << endl;
			dFactor->mbScore=move->getTargetMBScore();
			(*conditionSet)[cIter->first]=1;
			delete dFactor->potFunc;
			dFactor->potFunc=move->getDestPot();
			dFactor->updatePartialMeans(dFactor->potFunc->getAllPartialMeans());
			INTDBLMAP* plls=currPLLMap.begin()->second;
			(*plls)[dFactor->fId]=dFactor->mbScore;
			//cout << "### pot size after rm: " << regWts.size() << "," << edgeKey << endl;
			//if (tgtName == "YKL096W" || tgtName == "YGR249W" || tgtName == "YJL196C")
			//if (tgtName == "YKL096W")
			//{
			//	INTDBLMAP& regWts=dFactor->potFunc->getCondWeight();
			//	cout << "###\tRemoving " << tfName << " from " << tgtName << endl;
			//	cout << "###\tDeltaImprovement: " << move->getScoreImprovement() << endl;
			//	cout << "###\tCurrent edges:" << endl;
			//	for(INTDBLMAP_ITER mIter=regWts.begin();mIter!=regWts.end();mIter++)
			//	{
			//		string tname = varSet[mIter->first]->getName();
			//		double v = mIter->second;
			//		cout << "###\t" << tname << "\t" << tgtName << "\t" << v << endl;
			//	}
			//}
		}
	}
	//Get the module and update its indegree
	int mID=geneModuleID[v->getName()];
	map<string,int>* currIndegree=NULL;
	if(moduleIndegree.find(mID)==moduleIndegree.end())
	{
		currIndegree=new map<string,int>;
		moduleIndegree[mID]=currIndegree;
	}
	else
	{
		currIndegree=moduleIndegree[mID];
	}
	if(currIndegree->find(u->getName())==currIndegree->end())
	{
		//Shouldnt happen, because we are removing an existing edge!
	}
	else
	{	
		//cout <<"Updating regulator " << u->getName() <<" to module " << mID << endl;
		(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]-1;
	}
	if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
	{
		//Shouldnt happen, because we are removing an existing edge!
	}
	else
	{
		regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]-1;
	}
	//cout << "Made move for " << edgeKey.c_str() << " in condition " << csetind << endl;
	edgeRemoveUpdates[edgeKey]=1;
	(*conditionSet)[csetind]=1;
	int curriter=0;
	if(variableStatus.find(v->getName())==variableStatus.end())
	{
		variableStatus[v->getName()]=curriter;
	}
	else
	{
		variableStatus[v->getName()]=variableStatus[v->getName()]+1;
	}
	return 0;
}

int
GraphLearner::attemptMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
	VSET& varSet=varManager->getVariableSet();
	string edgeKey;
	Variable* u=varSet[move->getSrcVertex()];
	Variable* v=varSet[move->getTargetVertex()];
	edgeKey.append(u->getName().c_str());
	edgeKey.append("\t");
	edgeKey.append(v->getName().c_str());
	string tgtName = v->getName();
	string tfName  = u->getName();
	
	if(strcmp(u->getName().c_str(),"G114")==0 || strcmp(v->getName().c_str(),"G124")==0)
	{
		//cout <<"Stop here" << endl;
	}

	if(edgeConditionMap.find(edgeKey)==edgeConditionMap.end())
	{
		//cout <<"Edge " << edgeKey << " not found " << endl;
		return -1;
	}
	INTINTMAP* conditionSet=edgeConditionMap[edgeKey];
	int csetind=move->getConditionSetInd();
	if(affectedVars.find(1)==affectedVars.end())
	{
		cout << "No csvars for condition " <<csetind<< endl;
		exit(0);
	}
	INTINTMAP* csVars=affectedVars[1];
	if((csVars->find(move->getSrcVertex())!=csVars->end()) || (csVars->find(move->getTargetVertex())!=csVars->end()))
	{
		return -1;
	}
	addedEdges[edgeKey]=0;
	//(*csVars)[move->getSrcVertex()]=0;
	(*csVars)[move->getTargetVertex()]=0;
	INTINTMAP* condset=condsetMap[csetind];
	for(INTINTMAP_ITER cIter=condset->begin();cIter!=condset->end();cIter++)
	{
		if(cIter->second==0)
		{
			continue;
		}
		FactorGraph* csGraph=fgGraphSet[cIter->first];
		SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
		if(dFactor->mergedMB.find(move->getSrcVertex())!=dFactor->mergedMB.end())
		{
			//cout <<"Stop !! Trying to add the same edge " <<edgeKey << "   "<< v->getName() << endl;
		}
		//cout << "### Succeeded in add!" << endl;
		//cout << "### size before add: " << dFactor->mergedMB.size() << "," << edgeKey << endl;
		dFactor->mergedMB[move->getSrcVertex()]=0;
		//cout << "### size after add: " << dFactor->mergedMB.size() << "," << edgeKey << endl;
		dFactor->mbScore=move->getTargetMBScore();
		(*conditionSet)[cIter->first]=1;
		delete dFactor->potFunc;
		dFactor->potFunc=move->getDestPot();
		dFactor->updatePartialMeans(dFactor->potFunc->getAllPartialMeans());
		INTDBLMAP* plls=currPLLMap.begin()->second;
		(*plls)[dFactor->fId]=dFactor->mbScore;
		//if (tgtName == "YKL096W" || tgtName == "YGR249W" || tgtName == "YJL196C")
		//if (tgtName == "YKL096W")
		//{
		//	INTDBLMAP& regWts=dFactor->potFunc->getCondWeight();
		//	cout << "###\tAdding " << tfName << " to " << tgtName << endl;
		//	cout << "###\tDeltaImprovement: " << move->getScoreImprovement() << endl;
		//	cout << "###\tCurrent edges:" << endl;
		//	for(INTDBLMAP_ITER mIter=regWts.begin();mIter!=regWts.end();mIter++)
		//	{
		//		string tname = varSet[mIter->first]->getName();
		//		double v = mIter->second;
		//		cout << "###\t" << tname << "\t" << tgtName << "\t" << v << endl;
		//	}
		//}
	}
	//Get the module and update it's indegree
	int mID=geneModuleID[v->getName()];
	map<string,int>* currIndegree=NULL;
	if(moduleIndegree.find(mID)==moduleIndegree.end())
	{
		currIndegree=new map<string,int>;
		moduleIndegree[mID]=currIndegree;
	}
	else
	{
		currIndegree=moduleIndegree[mID];
	}
	if(currIndegree->find(u->getName())==currIndegree->end())
	{
		//cout <<"Adding new regulator " << u->getName() <<" to module " << mID << endl;
		(*currIndegree)[u->getName()]=1;
	}
	else
	{	
		//cout <<"Updating regulator " << u->getName() <<" to module " << mID << endl;
		(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]+1;
	}
	if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
	{
		regulatorModuleOutdegree[u->getName()]=1;
	}
	else
	{
		regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]+1;
	}
	//cout << "Made move for " << edgeKey.c_str() << " in condition " << csetind << endl;
	edgeUpdates[edgeKey]=1;
	(*conditionSet)[csetind]=1;
	int curriter=0;
	if(variableStatus.find(v->getName())==variableStatus.end())
	{
		variableStatus[v->getName()]=curriter;
	}
	else
	{
		variableStatus[v->getName()]=variableStatus[v->getName()]+1;
	}
	return 0;
}


int
GraphLearner::dumpAllGraphs(int currK,int foldid,int iter)
{
	VSET& varSet=varManager->getVariableSet();
	//char foldoutDirName[1024];
	char aFName[1024];
	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		string& dirname=outLocMap[eIter->first];
		//sprintf(aFName,"%s/var_mb_pw_k%d.txt",foldoutDirName,currK,iter);
		sprintf(aFName,"%s/prediction_k%d.txt",foldoutDirName,currK+1);
		ofstream oFile(aFName);
		//sprintf(aFName,"%s/bias.txt",foldoutDirName,currK);
		//ofstream bFile(aFName);
		for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
		{
			//Now we need to check if the subset of conditions corresponding to gIter->first
			//includes eIter->first
			INTINTMAP* cset=condsetMap[gIter->first];
			if((*cset)[eIter->first]==0)
			{
				continue;
			}
			FactorGraph* fg=gIter->second;
			fg->dumpVarMB_PairwiseFormat(oFile,varSet);
			//fg->dumpBias(bFile,varSet,reqdTargetList);
		}
		oFile.close();
	}
	for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
	{
		if(gIter->first>2)
		{
			char aFName[1024];
			string& dirname=outLocMap.begin()->second;
			//sprintf(aFName,"%s/fold%d/shared_mb_pw_k%d.txt",dirname.c_str(),foldid,currK);
			sprintf(aFName,"%s/fold%d/shared_mb_pw_k%d.txt",dirname.c_str(),foldid,currK+1);
			ofstream oFile(aFName);
			FactorGraph* fg=gIter->second;
			fg->dumpVarMB_PairwiseFormat(oFile,varSet);
			oFile.close();
		}
	}
	return 0;
}


int 
GraphLearner::getEdgeVars(string& edgeKey, char* var1, char* var2) 
{
	char buffer[1024];
	strcpy(buffer,edgeKey.c_str());
	char* pos=strchr(buffer,'\t');
	if(pos==NULL)
	{
		cout <<"Bad edge format" << endl;
		exit(0);
	}
	*pos='\0';
	int s=atoi(buffer);
	int d=atoi(pos+1);
	VSET& varSet=varManager->getVariableSet();
	strcpy(var1,varSet[s]->getName().c_str());
	strcpy(var2,varSet[d]->getName().c_str());
	return 0;
}


bool 
GraphLearner::checkMBSize(INTINTMAP* cset,int u,int v, int currK)
{
	bool check=true;
	for(INTINTMAP_ITER cIter=cset->begin();cIter!=cset->end();cIter++)
	{
		if(cIter->second==0)
		{
			continue;
		}
		FactorGraph* fg=fgGraphSet[cIter->first];
		SlimFactor* sFactor=fg->getFactorAt(u);
		SlimFactor* dFactor=fg->getFactorAt(v);
		if((dFactor->mergedMB.size()>=currK) && (dFactor->mergedMB.find(u)==dFactor->mergedMB.end()))
		{
			check=false;
		}
	}
	return check;
}


bool 
GraphLearner::checkMBSize(int cid, int u,int v, int currK)
{
	bool check=true;
	if(fgGraphSet.find(cid)==fgGraphSet.end())
	{
		return check;
	}
	FactorGraph* fg=fgGraphSet[cid];
	SlimFactor* sFactor=fg->getFactorAt(u);
	SlimFactor* dFactor=fg->getFactorAt(v);
	/*if((sFactor->mergedMB.size()>=currK) && (sFactor->mergedMB.find(v)==sFactor->mergedMB.end()))
	{
		check=false;
	}*/
	if((dFactor->mergedMB.size()>=currK) && (dFactor->mergedMB.find(u)==dFactor->mergedMB.end()))
	{
		check=false;
	}
	return check;
}

int
GraphLearner::populateGraphsFromInput()
{
	regulatorModuleOutdegree.clear();
	moduleIndegree.clear();
	VSET& varSet=varManager->getVariableSet();
	NET_T* inet = initnet->getNet();
	for(auto itr=inet->begin();itr!=inet->end();itr++)
	{
		string tfname = itr->first;
		map<string,double>* elist = itr->second;
		for (auto eitr=elist->begin(); eitr!=elist->end(); eitr++)
		{
			string tgname = eitr->first;
			Variable* u=NULL;
			Variable* v=NULL;
			INTINTMAP condSet;
			for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
			{
				condSet[eIter->first]=0;
			}

			int uid=varManager->getVarID(tfname);
			if (uid==-1)
			{
				continue;
			}
			u=varSet[uid];
			int vid=varManager->getVarID(tgname);
			if (vid==-1)
			{
				continue;
			}
			v=varSet[vid];
			if(geneModuleID.find(v->getName())==geneModuleID.end())
			{
				continue;
			}
			condSet[1]=1;

			int mID=geneModuleID[v->getName()];
			map<string,int>* currIndegree=NULL;
			if(moduleIndegree.find(mID)==moduleIndegree.end())
			{
				currIndegree=new map<string,int>;
				moduleIndegree[mID]=currIndegree;
			}
			else
			{
				currIndegree=moduleIndegree[mID];
			}
			if(currIndegree->find(u->getName())==currIndegree->end())
			{
				//cout <<"Adding new regulator " << u->getName() <<" to module " << mID << endl;
				(*currIndegree)[u->getName()]=1;
			}
			else
			{	
				//cout <<"Updating regulator " << u->getName() <<" to module " << mID << endl;
				(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]+1;
			}
			if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
			{
				regulatorModuleOutdegree[u->getName()]=1;
			}
			else
			{
				regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]+1;
			}

			for(INTINTMAP_ITER cIter=condSet.begin();cIter!=condSet.end();cIter++)
			{
				if(cIter->second==0)
				{
					continue;
				}
				FactorGraph* fg=fgGraphSet[cIter->first];
				SlimFactor* sFactor=fg->getFactorAt(u->getID());
				SlimFactor* dFactor=fg->getFactorAt(v->getID());
				//sFactor->mergedMB[dFactor->fId]=0;
				dFactor->mergedMB[sFactor->fId]=0;
				string edgeKey;
				edgeKey.append(u->getName().c_str());
				edgeKey.append("\t");
				edgeKey.append(v->getName().c_str());
				resumeEdges[edgeKey]=0;
			}
		}
	}
	//Populate all the potentials of all the graphs
	for(map<int,FactorGraph*>::iterator fIter=fgGraphSet.begin();fIter!=fgGraphSet.end();fIter++)
	{
		FactorGraph* condspecGraph=fIter->second;
		map<int,SlimFactor*>& factorSet=condspecGraph->getAllFactors();
		string& condKey=condsetIDKeyMap[fIter->first];
		PotentialManager* potMgr=potMgrSet[fIter->first];
		for(map<int,SlimFactor*>::iterator sIter=factorSet.begin();sIter!=factorSet.end();sIter++)
		{
			SlimFactor* sFactor=sIter->second;
			if (sFactor->potFunc != NULL)
			{
				delete sFactor->potFunc;
			}
			sFactor->potFunc=new Potential;
			sFactor->potFunc->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
			for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
			{
				sFactor->potFunc->setAssocVariable(varSet[mIter->first],Potential::MARKOV_BNKT);
				Variable* v = varSet[sFactor->fId];
				Variable* u = varSet[mIter->first];
				//cout << "populate!\t" << u->getName() << "\t" << v->getName() << endl;
			}
			sFactor->potFunc->potZeroInit();
			potMgr->populatePotential(sFactor->potFunc,random);
			sFactor->potFunc->initMBCovMean();
		}
	}
	for(map<int,FactorGraph*>::iterator fIter=fgGraphSet.begin();fIter!=fgGraphSet.end();fIter++)
	{
		FactorGraph* condspecGraph=fIter->second;
		map<int,SlimFactor*>& factorSet=condspecGraph->getAllFactors();
		string& condKey=condsetIDKeyMap[fIter->first];
		for(map<int,SlimFactor*>::iterator sIter=factorSet.begin();sIter!=factorSet.end();sIter++)
		{
			SlimFactor* sFactor=sIter->second;
			string gname = varSet[sFactor->vIds[0]]->getName();
			//if (gname == "YKL096W" || gname == "YGR249W" || gname == "YJL196C")
			//if (gname == "YKL096W")
			//{
			//	INTDBLMAP& regWts=sFactor->potFunc->getCondWeight();
			//	for(INTDBLMAP_ITER mIter=regWts.begin();mIter!=regWts.end();mIter++)
			//	{
			//		string tname = varSet[mIter->first]->getName();
			//		double v = mIter->second;
			//		cout << "###\tINITIALIZED\t" << tname << "\t" << gname << "\t" << v << endl;
			//	}
			//}
		}
	}
	return 0;
}

int
GraphLearner::populateGraphsFromFile(const char* aFName)
{
	regulatorModuleOutdegree.clear();
	moduleIndegree.clear();
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=varManager->getVariableSet();
	while(inFile.good())
	{	
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Variable* u=NULL;
		Variable* v=NULL;
		INTINTMAP condSet;
		for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
		{
			condSet[eIter->first]=1;
		}
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				int vid=varManager->getVarID(tok);
				u=varSet[vid];
			}
			else if(tokCnt==1)
			{
				int vid=varManager->getVarID(tok);
				v=varSet[vid];
			}
			/*
			else
			{
				int cid=atoi(tok);
				condSet[cid]=1;
			}
			*/
			tok=strtok(NULL,"\t");
			tokCnt++;
		}

		int mID=geneModuleID[v->getName()];
		map<string,int>* currIndegree=NULL;
		if(moduleIndegree.find(mID)==moduleIndegree.end())
		{
			currIndegree=new map<string,int>;
			moduleIndegree[mID]=currIndegree;
		}
		else
		{
			currIndegree=moduleIndegree[mID];
		}
		if(currIndegree->find(u->getName())==currIndegree->end())
		{
			//cout <<"Adding new regulator " << u->getName() <<" to module " << mID << endl;
			(*currIndegree)[u->getName()]=1;
		}
		else
		{	
			//cout <<"Updating regulator " << u->getName() <<" to module " << mID << endl;
			(*currIndegree)[u->getName()]=(*currIndegree)[u->getName()]+1;
		}
		if(regulatorModuleOutdegree.find(u->getName())==regulatorModuleOutdegree.end())
		{
			regulatorModuleOutdegree[u->getName()]=1;
		}
		else
		{
			regulatorModuleOutdegree[u->getName()]=regulatorModuleOutdegree[u->getName()]+1;
		}

		/*string condKey;
		genCondSetKey(condSet,condKey);
		if(pooledData.find(condKey)==pooledData.end())
		{
			genPooledDataset(condSet,condKey);
		}
		PotentialManager* potMgr=pooledPotentials[condKey];
		int csetid=condsetKeyIDMap[condKey];*/
		for(INTINTMAP_ITER cIter=condSet.begin();cIter!=condSet.end();cIter++)
		{
			if(cIter->second==0)
			{
				continue;
			}
			FactorGraph* fg=fgGraphSet[cIter->first];
			SlimFactor* sFactor=fg->getFactorAt(u->getID());
			SlimFactor* dFactor=fg->getFactorAt(v->getID());
			//sFactor->mergedMB[dFactor->fId]=0;
			dFactor->mergedMB[sFactor->fId]=0;
			string edgeKey;
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());
			resumeEdges[edgeKey]=0;
		}
	}
	//Populate all the potentials of all the graphs
	for(map<int,FactorGraph*>::iterator fIter=fgGraphSet.begin();fIter!=fgGraphSet.end();fIter++)
	{
		FactorGraph* condspecGraph=fIter->second;
		map<int,SlimFactor*>& factorSet=condspecGraph->getAllFactors();
		string& condKey=condsetIDKeyMap[fIter->first];
		PotentialManager* potMgr=potMgrSet[fIter->first];
		for(map<int,SlimFactor*>::iterator sIter=factorSet.begin();sIter!=factorSet.end();sIter++)
		{
			SlimFactor* sFactor=sIter->second;
			if (sFactor->potFunc != NULL)
			{
				delete sFactor->potFunc;
			}
			sFactor->potFunc=new Potential;
			sFactor->potFunc->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
			for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
			{
				sFactor->potFunc->setAssocVariable(varSet[mIter->first],Potential::MARKOV_BNKT);
			}
			sFactor->potFunc->potZeroInit();
			potMgr->populatePotential(sFactor->potFunc,random);
			sFactor->potFunc->initMBCovMean();
		}
	}
	inFile.close();
	return 0;
}

int
GraphLearner::normalizeWeight(INTDBLMAP& wtVect)
{
	double maxExp=-1000000000;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		if(dIter->second>maxExp)
		{
			maxExp=dIter->second;
		}
	}
	double pow=-1*maxExp;
	double total=0;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		total=total+exp(dIter->second+pow);
	}
	total=log(total)-pow;
	//correction factor
	double newtotal=0;
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		double scWt=1e-5+exp(tempval-total);
		dIter->second=scWt;
		newtotal=newtotal+scWt;
	}
	for(INTDBLMAP_ITER dIter=wtVect.begin();dIter!=wtVect.end();dIter++)
	{
		double tempval=dIter->second;
		double scWt=tempval/newtotal;
		dIter->second=scWt;
	}
	return 0;
}


int
GraphLearner::updatePotentials()
{
	VSET& varSet=varManager->getVariableSet();
	/*for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
	{
		FactorGraph* graph=gIter->second;
		map<int,SlimFactor*>& factorSet=graph->getAllFactors();
		string& condKey=condsetIDKeyMap[gIter->first];
		PotentialManager* potMgr=pooledPotentials[condKey];
		for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
		{
			SlimFactor* sFactor=fIter->second;
			potMgr->populatePotential(sFactor->potFunc,false);
			sFactor->potFunc->initMBCovMean();
			if(sFactor->potFunc->getCondVariance()<=0)
			{
				cout <<"Found negative covariance for "<< sFactor->fId << " "  << varSet[sFactor->fId]->getName() <<  " in condition " << gIter->first  << endl;
				sFactor->potFunc->setCondVariance(1e-20);
			}
		}
	}*/

	//Now reestimate the plls for each dataset and potential
	INTDBLMAP* plls=currPLLMap[evMgrSet.begin()->first];
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* var=vIter->second;
		double newPLL_s=getNewPLLScore_Condition(-1,vIter->first,NULL);
		(*plls)[vIter->first]=newPLL_s;
	}

	return 0;
}

int
GraphLearner::initPhysicalDegree()
{
	for(map<int,map<string,int>*>::iterator mIter=moduleGeneSet.begin();mIter!=moduleGeneSet.end();mIter++)
	{
		map<string,int>* indegree=NULL;
		map<string,map<string,int>*> innet;
		map<string,int>* geneSet=mIter->second;
		for(map<string,map<string,map<string,double>*>*>::iterator gpIter=priorgraphmap.begin();gpIter!=priorgraphmap.end();gpIter++)
		{
			map<string,int>enrichedTFs;
			map<string,map<string,double>*>* priorgraph = gpIter->second;
			getEnrichedTFs(enrichedTFs,geneSet,*priorgraph);
			for(map<string,int>::iterator tfIter=enrichedTFs.begin();tfIter!=enrichedTFs.end();tfIter++)
			{
				map<string,double>* motiftgts=(*priorgraph)[tfIter->first];
				map<string,int>* ttgts;
				if (innet.find(tfIter->first) == innet.end())
				{
					ttgts = new map<string,int>();
					innet[tfIter->first] = ttgts;
				}
				else
				{
					ttgts = innet[tfIter->first];
				}
				for(map<string,double>::iterator gIter=motiftgts->begin();gIter!=motiftgts->end();gIter++)
				{
					if(geneSet->find(gIter->first)==geneSet->end())
					{
						continue;
					}
					(*ttgts)[gIter->first] = 0;
				}
			}
		}
		for (map<string,map<string,int>*>::iterator tItr=innet.begin();tItr!=innet.end();tItr++)
		{
			string tf = tItr->first;
			map<string,int>* ttgts = tItr->second;
			if (ttgts->size() == 0)
			{
				continue;
			}
			if (indegree == NULL)
			{
				indegree=new map<string,int>;
			}
			(*indegree)[tf] = ttgts->size();
		}
		
		if(indegree!=NULL)
		{
			moduleIndegree[mIter->first]=indegree;
			//cout <<"Physical_indegree for module " << mIter->first << endl;
			for(map<string,int>::iterator dIter=indegree->begin();dIter!=indegree->end();dIter++)
			{
				//cout << dIter->first <<"\t" << dIter->second << endl;
				if(regulatorModuleOutdegree.find(dIter->first)==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[dIter->first]=dIter->second;
				}
				else
				{
					regulatorModuleOutdegree[dIter->first]=regulatorModuleOutdegree[dIter->first]+dIter->second;
				}
			}
		}
	}
	return 0;
}

int
GraphLearner::getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet)
{
	VSET& varSet=varManager->getVariableSet();
	int total=varSet.size();
	int k=genes->size();
	HyperGeomPval hgp;
	for(map<string,map<string,double>*>::iterator fIter=edgeSet.begin();fIter!=edgeSet.end();fIter++)
	{
		int vID=varManager->getVarID(fIter->first.c_str());
		if(vID<0)
		{
			continue;
		}
		map<string,double>* tgtSet=fIter->second;
		//int n=tgtSet->size();
		int n=0;
		int hit=0;
		for(map<string,double>::iterator gIter=tgtSet->begin();gIter!=tgtSet->end();gIter++)
		//for(map<string,int>::iterator gIter=genes->begin();gIter!=genes->end();gIter++)
		{
			
			int vID=varManager->getVarID(gIter->first.c_str());
                        if(vID<0)
                        {
                                continue;
                        }
                        n++;
			//if(tgtSet->find(gIter->first)==tgtSet->end())
                	if(genes->find(gIter->first)==genes->end())
			{
				continue;
			}
			hit++;
		}
		double enpval=hgp.getOverRepPval(k,hit,n,total-n);
		if(enpval<0.05 && hit>4)
		//if(hit>0)
		{
			tfSet[fIter->first]=hit;
		}
	}
	return 0;
}

double
GraphLearner::getModuleContrib(string& tgtName, string& tfName)
{
	//return 0;
	int moduleID=-1;
	if(geneModuleID.find(tgtName)==geneModuleID.end())
	{
		return log(1.0/(double)geneModuleID.size());	
	}
	moduleID=geneModuleID[tgtName];
	int moduleSize=moduleGeneSet[moduleID]->size();
	double contrib=0;
	if(moduleIndegree.find(moduleID)==moduleIndegree.end())
	{
		contrib=1.0/(moduleSize);
		return log(contrib);
	}
	int degree=0;
	double total=0;
	map<string,int>* moddegree=moduleIndegree[moduleID];
	for(map<string,int>::iterator rIter=moddegree->begin();rIter!=moddegree->end();rIter++)
	{
		total=total+rIter->second;
	}
	if(moddegree->find(tfName)==moddegree->end())
	{
		contrib=1.0/(moduleSize);
		return log(contrib);		
	}
	degree=(*moddegree)[tfName];
	contrib=log((double)(degree)/(moduleSize));
	return contrib;
}


double
GraphLearner::getModuleContribLogistic(string& tgtName, string& tfName)
{
	double mbeta2=beta_motif;
	if(geneModuleID.find(tgtName)==geneModuleID.end())
	{
		return 0;
	}

	int regDegree=0;
	if(regulatorModuleOutdegree.find(tfName)!=regulatorModuleOutdegree.end())
	{
		regDegree=regulatorModuleOutdegree[tfName];
	}
	int moduleID=geneModuleID[tgtName];
	if(moduleIndegree.find(moduleID)==moduleIndegree.end())
	{
		return 0;
	}
	int degree=0;
	double total=0;
	map<string,int>* moddegree=moduleIndegree[moduleID];
	if(moddegree->find(tfName)==moddegree->end())
	{
		return 0;
	}
	degree=(*moddegree)[tfName];
	double contrib=((double) degree)/((double) regDegree);
	double mprior2=mbeta2*contrib;
	return mprior2;
}

//To redefine the modules we will start with the original set of modules 
//For each original module, find for every gene its pairwise similarity to every other
//gene. merge two nodes that have the greatest pairwise similarity. replace by the merged
//regulatory program. recompute similarity of all nodes to this merged node. repeat with
//finding the next most similar pair of nodes.

//To redefine the modules we will start with the original set of modules 
//For each original module, find for every gene its pairwise similarity to every other
//gene. merge two nodes that have the greatest pairwise similarity. replace by the merged
//regulatory program. recompute similarity of all nodes to this merged node. repeat with
//finding the next most similar pair of nodes.
int
GraphLearner::redefineModules()
{
	FactorGraph* aGraph=fgGraphSet.begin()->second;
	map<int,map<string,int>*> newModules;
	for(map<int,map<string,int>*>::iterator gIter=moduleGeneSet.begin();gIter!=moduleGeneSet.end();gIter++)
	{
		map<string,int>* moduleMembers=gIter->second;
		map<string,HierarchicalClusterNode*> nodeSet;
		for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
		{
			int mID=varManager->getVarID(mIter->first.c_str());
			SlimFactor* mFactor=aGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;
			HierarchicalClusterNode* node=new HierarchicalClusterNode;
			for(INTINTMAP_ITER bIter=mbvars1.begin();bIter!=mbvars1.end();bIter++)
			{
				node->attrib[bIter->first]=bIter->second;
			}
			node->nodeName.append(mIter->first);
			nodeSet[mIter->first]=node;
		}
		HierarchicalCluster hc;
		hc.cluster(newModules,nodeSet,clusterThreshold);
	}
	moduleGeneSet.clear();
	geneModuleID.clear();
	for(map<int,map<string,int>*>::iterator mIter=moduleIndegree.begin();mIter!=moduleIndegree.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	moduleIndegree.clear();
	regulatorModuleOutdegree.clear();
	VSET& varSet=varManager->getVariableSet();
	for(map<int,map<string,int>*>::iterator mIter=newModules.begin();mIter!=newModules.end();mIter++)
	{
		moduleGeneSet[mIter->first]=mIter->second;
		map<string,int>* geneSet=mIter->second;
		map<string,int>* indegree=new map<string,int>;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			geneModuleID[gIter->first]=mIter->first;
			int mID=varManager->getVarID(gIter->first.c_str());
			SlimFactor* mFactor=aGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;
		
			for(INTINTMAP_ITER nIter=mbvars1.begin();nIter!=mbvars1.end();nIter++)
			{
				Variable* var=varSet[nIter->first];
				if(indegree->find(var->getName())==indegree->end())
				{
					(*indegree)[var->getName()]=1;
				}
				else
				{
					(*indegree)[var->getName()]=(*indegree)[var->getName()]+1;
				}
				if(regulatorModuleOutdegree.find(var->getName())==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[var->getName()]=1;
				}	
				else
				{
					regulatorModuleOutdegree[var->getName()]=regulatorModuleOutdegree[var->getName()]+1;
				}
			}
		}
		moduleIndegree[mIter->first]=indegree;
	}

	return 0;
}

int
GraphLearner::redefineModules_Global()
{
	FactorGraph* aGraph=fgGraphSet.begin()->second;
	map<int,map<string,int>*> newModules;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,int> genesWithNoNeighbors;
	EvidenceManager* evMgr=evMgrSet.begin()->second;
	INTINTMAP& tSet=evMgr->getTrainingSet();
	for(map<int,map<string,int>*>::iterator gIter=moduleGeneSet.begin();gIter!=moduleGeneSet.end();gIter++)
	{
		map<string,int>* moduleMembers=gIter->second;
		for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
		{
			int mID=varManager->getVarID(mIter->first.c_str());
			if(mID<0)
			{
				continue;
			}
			SlimFactor* mFactor=aGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;
			INTDBLMAP& regWts=mFactor->potFunc->getCondWeight();
			if(mbvars1.size()==0)
			{
				genesWithNoNeighbors[mIter->first]=0;
				continue;
			}

			HierarchicalClusterNode* node=new HierarchicalClusterNode;
			//for(INTINTMAP_ITER bIter=mbvars1.begin();bIter!=mbvars1.end();bIter++)
			for(INTDBLMAP_ITER bIter=regWts.begin();bIter!=regWts.end();bIter++)
			{
				node->attrib[bIter->first]=bIter->second;
			}
			//Get the expression data
			for(INTINTMAP_ITER eIter=tSet.begin();eIter!=tSet.end();eIter++)
			{	
				EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
				Evidence* evid=(*evidMap)[mID];
				double v=evid->getEvidVal();
				node->expr.push_back(v);
			}
			node->nodeName.append(mIter->first);
			nodeSet[mIter->first]=node;
			node->size=1;
		}
	}
	HierarchicalCluster hc;
	hc.setOutputDir(foldoutDirName);
	hc.setVariableManager(varManager);
	hc.cluster(newModules,nodeSet,clusterThreshold);
	moduleGeneSet.clear();
	geneModuleID.clear();
	regulatorModuleOutdegree.clear();
	for(map<int,map<string,int>*>::iterator mIter=moduleIndegree.begin();mIter!=moduleIndegree.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	moduleIndegree.clear();
	VSET& varSet=varManager->getVariableSet();
	int largestModuleID=0;
	char moduleFName[1024];
	
	string& dirname=outLocMap[evMgrSet.begin()->first];
	sprintf(moduleFName,"%s/fold%d/modules.txt",dirname.c_str(),currFold);
	ofstream modFile(moduleFName);
	for(map<int,map<string,int>*>::iterator mIter=newModules.begin();mIter!=newModules.end();mIter++)
	{
		moduleGeneSet[mIter->first]=mIter->second;
		map<string,int>* geneSet=mIter->second;
		map<string,int>* indegree=new map<string,int>;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			modFile << gIter->first <<"\t" << mIter->first << endl;
			geneModuleID[gIter->first]=mIter->first;
			int mID=varManager->getVarID(gIter->first.c_str());
			SlimFactor* mFactor=aGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;
		
			for(INTINTMAP_ITER nIter=mbvars1.begin();nIter!=mbvars1.end();nIter++)
			{
				Variable* var=varSet[nIter->first];
				if(indegree->find(var->getName())==indegree->end())
				{
					(*indegree)[var->getName()]=1;
				}
				else
				{
					(*indegree)[var->getName()]=(*indegree)[var->getName()]+1;
				}
				if(regulatorModuleOutdegree.find(var->getName())==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[var->getName()]=1;
				}	
				else
				{
					regulatorModuleOutdegree[var->getName()]=regulatorModuleOutdegree[var->getName()]+1;
				}
			}
		}
		moduleIndegree[mIter->first]=indegree;
		largestModuleID=mIter->first;
	}
	modFile.close();
	for(map<string,int>::iterator gIter=genesWithNoNeighbors.begin();gIter!=genesWithNoNeighbors.end();gIter++)
	{
		largestModuleID++;
		map<string,int>* newmodule=new map<string,int>;
		(*newmodule)[gIter->first]=0;
		moduleGeneSet[largestModuleID]=newmodule;
		geneModuleID[gIter->first]=largestModuleID;
	}
	genesWithNoNeighbors.clear();
	return 0;
}

int
GraphLearner::setRestrictedList(vector<string>* names)
{
	restrictedVarList.clear();
	for (int i=0; i<names->size(); i++)
	{
		restrictedVarList[names->at(i)] = 0;
	}
	return 0;
}

int
GraphLearner::setRestrictedList(VariableManager* tempMgr)
{
	restrictedVarList.clear();
	VSET vset = tempMgr->getVariableSet();
	for (auto itr=vset.begin(); itr!=vset.end(); itr++)
	{
		Variable* var = itr->second;
		restrictedVarList[var->getName()] = 0;
	}
	return 0;
}

int
GraphLearner::setDefaultModuleMembership(map<string,int>& moduleID)
{
	geneModuleID.clear();
	moduleGeneSet.clear();
	for (auto itr=moduleID.begin();itr!=moduleID.end();itr++)
	{
		string g = itr->first;
		int v = itr->second;
		geneModuleID[g] = v;
		map<string,int>* geneSet=NULL;
		if (moduleGeneSet.find(v)==moduleGeneSet.end())
		{
			geneSet=new map<string,int>;
			moduleGeneSet[v]=geneSet;
		}
		else
		{
			geneSet = moduleGeneSet[v];
		}
		(*geneSet)[g] = 0;
	}
	return 0;
}

int
GraphLearner::setDefaultModuleMembership(vector<string>* names)
{
	map<string,int>* geneSet=NULL;
	geneSet=new map<string,int>;
	moduleGeneSet[0]=geneSet;

	for(int itr=0;itr<names->size();itr++)
	{
		(*geneSet)[names->at(itr)]=0;
		geneModuleID[names->at(itr)]=0;
	}
	return 0;
}

int
GraphLearner::setDefaultModuleMembership(VariableManager* tempVarMngr)
{
	VSET vset = tempVarMngr->getVariableSet();
	map<string,int>* geneSet=NULL;
	geneSet=new map<string,int>;
	moduleGeneSet[0]=geneSet;

	for(auto itr=vset.begin();itr!=vset.end();itr++)
	{
		Variable* v=itr->second;
		(*geneSet)[v->getName()]=0;
		geneModuleID[v->getName()]=0;
	}
	return 0;
}

int
GraphLearner::setPriorGraph_All(const char* netname, EdgeList* pnet, double b)
{
	betamap[netname] = b;
	priorgraphmap[netname] = pnet->getNet();
	priorpointers[netname] = pnet;
	return 0;
}

EdgeList*
GraphLearner::getLearnedGraph()
{
	EdgeList* onet = new EdgeList;
	VSET& varSet=varManager->getVariableSet();
	for(map<int,EvidenceManager*>::iterator eIter=evMgrSet.begin();eIter!=evMgrSet.end();eIter++)
	{
		for(map<int,FactorGraph*>::iterator gIter=fgGraphSet.begin();gIter!=fgGraphSet.end();gIter++)
		{
			//Now we need to check if the subset of conditions corresponding to gIter->first
			//includes eIter->first
			INTINTMAP* cset=condsetMap[gIter->first];
			if((*cset)[eIter->first]==0)
			{
				continue;
			}
			FactorGraph* fg=gIter->second;
			fg->getGraph(onet,varSet);
		}
	}
	return onet;
}

map<string,int>* 
GraphLearner::getLearnedModule()
{
	map<string,int>* omod = new map<string,int>;
	for (auto itr=geneModuleID.begin();itr!=geneModuleID.end();itr++)
	{
		(*omod)[itr->first] = itr->second;
	}
	return omod;
}

int 
GraphLearner::writeIterCnt(int currFold, int iter)
{
	string& dirname=outLocMap[evMgrSet.begin()->first];
	char aFName[1024];
	sprintf(aFName,"%s/fold%d/iter.txt",dirname.c_str(),currFold);
	ofstream outFile(aFName);
	outFile << iter << endl;
	outFile.close();
	return 0;
}

int
GraphLearner::readIterCnt(const char* aFName)
{
	int iter = 0;
	ifstream inFile(aFName);
	inFile >> iter;
	inFile.close();
	return iter;
}

int
GraphLearner::doTar()
{
	char tarCmd[1024];
	sprintf(tarCmd,"tar cvzf %s.tar.gz %s",topOutputDirName,topOutputDirName);
	system(tarCmd);
	return 0;
}

int
GraphLearner::writeAddRemEdges(const char* aFName, map<string,int>& emap)
{
	ofstream outFile(aFName);
	for (map<string,int>::iterator eitr=emap.begin();eitr!=emap.end();eitr++)
	{
		outFile << eitr->first << endl;
	}
	outFile.close();
	return 0;
}

int
GraphLearner::readAddRemEdges(const char* aFName, map<string,int>& emap)
{
	ifstream inFile(aFName);
	char buffer[1024];
	VSET& varSet=varManager->getVariableSet();
	while(inFile.good())
	{	
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Variable* u=NULL;
		Variable* v=NULL;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				int vid=varManager->getVarID(tok);
				u=varSet[vid];
			}
			else if(tokCnt==1)
			{
				int vid=varManager->getVarID(tok);
				v=varSet[vid];
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		string edgeKey;
		edgeKey.append(u->getName().c_str());
		edgeKey.append("\t");
		edgeKey.append(v->getName().c_str());
		emap[edgeKey]=0;
	}
	inFile.close();
	return 0;
}

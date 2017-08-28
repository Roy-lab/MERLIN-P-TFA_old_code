#include "Framework.H"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <pthread.h>
#include <semaphore.h>

#include "ThreadManager.H"
#include "MemoryCheck.H"

Framework::Framework()
{
	regMngr = NULL;
	tgtMngr = NULL;
	evidMngr = NULL;
	priorNet = NULL;
	lambda = 0;
	CV=0;
	sparsity=-5;
	beta=5;
	numNcaRep=10;
	numThread=1;
	lShouldResume = false;
	lResumeMERLIN = false;
	resIter = -1;
}

Framework::~Framework()
{
	if (regMngr!=NULL)
		delete regMngr;
	if (tgtMngr!=NULL)
		delete tgtMngr;
	if (evidMngr!=NULL)
		delete evidMngr;
	if (priorNet!=NULL)
		delete priorNet;
}

bool
Framework::shouldResume()
{
	return lShouldResume;
}

GraphLearner*
Framework::getGL(int itrNum, int subNum, map<string,int>& geneModuleID, EdgeList* initNet, bool sl)
{
	GraphLearner* glearner = new GraphLearner;
	
	VariableManager* vm = allGenesMngr->copyMe();
	Matrix* dmat = allEvidMngr->getDataMat();

	EvidenceManager* emgr = new EvidenceManager;
	emgr->setVariableManager(vm);
	emgr->setDataMat(dmat);
	delete dmat;
	glearner->setVariableManager(vm);
	glearner->setGlobalEvidenceManager(emgr);

	char od[1025];
	sprintf(od,"%s/%d/%d/",outdir,itrNum,subNum);
	glearner->setOutputDirName(od);	
	glearner->setMaxFactorSize(1);
	glearner->setMaxFactorSize_Approx(300);

	glearner->setRestrictedList(regMngr);

	//glearner->setDefaultModuleMembership(tgtMngr);
	glearner->setDefaultModuleMembership(geneModuleID);

	glearner->setClusteringThreshold(0.6);
	//glearner->setBeta1(-5);
	glearner->setBeta1(sparsity);
	//glearner->setBeta_Motif(4);
	glearner->setBeta_Motif(4);//This is modularity

	glearner->setShouldLoad(sl);
	glearner->setTopOutputDirName(outdir);

	EdgeList* cnet = pnet->copyMe();
	//glearner->setPriorGraph_All("motif",cnet,5);
	glearner->setPriorGraph_All("motif",cnet,beta);
	if (initNet!=NULL)
	{
		EdgeList* tnet = initNet->copyMe();
		glearner->setInitNet(tnet);
	}

	glearner->initPartitions(1);

	return glearner;
}

void*
Framework::runOneGL(void* ptr)
{
	void** pptr = (void**)ptr;
	ThreadManager* tm = (ThreadManager*) pptr[0];

	//fprintf(stdout,"startRunOne\n");
	GraphLearner* glearner = (GraphLearner*) pptr[1];
	glearner->doCrossValidation(1);
	EdgeList* onet = glearner->getLearnedGraph();
	map<string,int>* omod = glearner->getLearnedModule();
	delete glearner;

	pthread_mutex_t& qlock = tm->getLock();
	sem_t& qsem = tm->getSem();

	pthread_mutex_lock(&qlock);
	vector<EdgeList*>* outnets = (vector<EdgeList*>*) pptr[2];
	vector<map<string,int>* >* outmods = (vector<map<string,int>* >*) pptr[3];
	outnets->push_back(onet);
	outmods->push_back(omod);
	pthread_mutex_unlock(&qlock);
	sem_post(&qsem);
	return NULL;
}

EdgeList*
Framework::inferNetwork(int itrNum, EdgeList* initNet)
{
	ThreadManager* thMgr = new ThreadManager(numThread,&runOneGL);
	VSET vset = tgtMngr->getVariableSet();
	int vsize = vset.size();
	//If I have more threads than genes, I will use less threads
	if (numThread>vsize)
	{
		numThread = vsize;
	}
	int rsize = (int)ceil(double(vsize)/double(numThread));
	vector<EdgeList*>* outnets = new vector<EdgeList*>;
	vector<map<string,int>* >* outmods = new vector<map<string,int>* >;
	outnets->clear();
	for (int ritr=0;ritr<vsize;ritr+=rsize)
	{
		int b=ritr;
		int e=min(ritr+rsize,vsize);
		map<string,int> names;
		names.clear();
		for (int j=b;j<e;j++)
		{
			names[vset[j]->getName()] = moduleID[vset[j]->getName()];
		}
		GraphLearner* glearner = getGL(itrNum,ritr/rsize,names,initNet,lResumeMERLIN);
		void** pptr = new void*[4];
		pptr[0] = (void*) thMgr;
		pptr[1] = (void*) glearner;
		pptr[2] = (void*) outnets;
		pptr[3] = (void*) outmods;
		thMgr->addInput((void*)pptr);
	}
	thMgr->run();
	delete thMgr;
	EdgeList* unewNet = new EdgeList;
	unewNet->setUnion(outnets);
	EdgeList* newNet = unewNet->percentileRank();
	delete unewNet;
	moduleID.clear();
	for (int i=0;i<outmods->size();i++)
	{
		map<string,int>* omod = outmods->at(i);
		for (auto mitr=omod->begin();mitr!=omod->end();mitr++)
		{
			moduleID[mitr->first] = mitr->second;
		}
		delete omod;
	}
	outmods->clear();
	char modoname[1024];
	sprintf(modoname,"%s/module_%d.txt",outdir,itrNum);
	writeClusters(modoname,moduleID);

	char netoname[1024];
	sprintf(netoname,"%s/merged_%d.txt",outdir,itrNum);
	newNet->writeNet(netoname);
	for (int i=0;i<outnets->size();i++)
	{
		EdgeList* e = outnets->at(i);
		delete e;
	}
	delete outnets;

	//TAR!!!
	char tarCmd[1024];
	sprintf(tarCmd,"tar cvzf %s.tar.gz %s",outdir,outdir);
	system(tarCmd);

	return newNet;
}

int
Framework::inferNCA(int itrNum, EdgeList* initNet)
{
	//In case of duplicate (both expression and TFA regulator),
	//I am taking Max. Maybe I should take mean?
	EdgeList* initNet2 = initNet->removeSuffix("_nca");
	char ioutdir[1025];
	sprintf(ioutdir,"%s/%d/",outdir,itrNum);
	Graph* updatedPrior = new Graph;
	updatedPrior->setVariableManagers(regMngr,tgtMngr);
	updatedPrior->setNet(initNet2);
	//NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, priorNet, lambda, 10, numThread, ioutdir);
	NCALearner* ncalrnr = new NCALearner(regMngr, tgtMngr, evidMngr, updatedPrior, lambda, CV, numNcaRep, numThread, ioutdir);
	ncalrnr->start();
	VariableManager* nregMngr;
	EvidenceManager* ntfaMngr;
	ncalrnr->getData(nregMngr,ntfaMngr,"_nca");
	delete ncalrnr;

	regMngr->mergeVarSets(nregMngr);
	allGenesMngr->mergeVarSets(nregMngr);
	allEvidMngr->updateEvidence(ntfaMngr);

	delete ntfaMngr;
	delete initNet2;

	return 0;
}


bool
Framework::isConverged(int itrNum, EdgeList* oldNet, EdgeList* newNet)
{
	int oldcount = oldNet->countEdges();
	int newcount = newNet->countEdges();
	int overlap  = newNet->getOverlap(oldNet);
	double jacc = double(overlap)/double(oldcount+newcount-overlap);
	cout << "Iteration: " << itrNum << ", old edge count: " << oldcount << ", new edge count: " << newcount;
	cout << ", overlap: " << overlap <<  ", jaccard: " << jacc << endl;
	if (itrNum>4)// || jacc>0.99)
	{
		return true;
	}
	return false;
}

bool
Framework::doNextStep(int itrNum, EdgeList*& newNet)
{
	MemoryCheck mc;
	char message[1024];
	EdgeList* oldNet = NULL;
	oldNet = newNet;
	mc.begin();
	inferNCA(itrNum,oldNet);
	mc.end();
	sprintf(message,"inferNCA(%d)",itrNum);
	mc.print(message);

	mc.begin();
	newNet = inferNetwork(itrNum,oldNet);
	mc.end();
	sprintf(message,"inferNetwork(%d)",itrNum);
	mc.print(message);

	bool conv = isConverged(itrNum,oldNet,newNet);
	
	delete oldNet;
	return conv;
}

int
Framework::start()
{
	//Infer one network
	//while True
	//{
	//	Infer NCA
	//	Infer one network
	//}

	int itrNum=0;

	MemoryCheck mc;
	//Start with MotifTFA
	mc.begin();
	inferNCA(itrNum,pnet);
	mc.end();
	mc.print("inferNCA(0)");

	mc.begin();
	EdgeList* newNet = NULL;
	newNet = inferNetwork(itrNum, NULL);
	EdgeList* oldNet = NULL;
	mc.end();
	mc.print("inferNetwork(0)");

	//return 0;
	bool conv=false;

	while(true)
	{
		itrNum++;
		conv = doNextStep(itrNum,newNet);
		if (conv)
		{
			break;
		}
	}
	delete newNet;
	return 0;
}

int
Framework::readTFA(int itrNum)
{
	char tfaname[1024];
	sprintf(tfaname,"%s/%d/tfa.txt",outdir,itrNum);
	VariableManager* oregMngr = new VariableManager;
	oregMngr->readVariables(tfaname);
	EvidenceManager* otfaMngr = new EvidenceManager;
	otfaMngr->setVariableManager(oregMngr);
	otfaMngr->loadEvidenceFromFile_Continuous(tfaname);

	VariableManager* nregMngr = oregMngr->copyMeWithSuffix("_nca");
	EvidenceManager* ntfaMngr = new EvidenceManager;
	ntfaMngr->setVariableManager(nregMngr);
	Matrix* mat = otfaMngr->getDataMat();
	ntfaMngr->setDataMat(mat);
	delete mat;

	regMngr->mergeVarSets(nregMngr);
	allGenesMngr->mergeVarSets(nregMngr);
	allEvidMngr->updateEvidence(ntfaMngr);
	delete ntfaMngr;
	delete otfaMngr;
	return 0;
}

int
Framework::resume()
{
	int itrNum=resIter;
	EdgeList* newNet = new EdgeList;
	if(lResumeMERLIN)
	{
		readTFA(itrNum);
		EdgeList* oldNet = NULL;
		if(itrNum>0)
		{
			oldNet = new EdgeList;
			char netoname[1024];
			sprintf(netoname,"%s/merged_%d.txt",outdir,itrNum-1);
			oldNet->readNet(netoname);
		}
		newNet = inferNetwork(itrNum,oldNet);

		if (oldNet!=NULL)
		{
			bool conv = isConverged(itrNum,oldNet,newNet);
			if(conv)
			{
				delete newNet;
				delete oldNet;
				return 0;
			}
		}
		lResumeMERLIN=false;
	}
	else
	{
		char modoname[1024];
		sprintf(modoname,"%s/module_%d.txt",outdir,itrNum);
		readClusters(modoname,moduleID);

		char netoname[1024];
		sprintf(netoname,"%s/merged_%d.txt",outdir,itrNum);
		newNet->readNet(netoname);
	}

	bool conv=false;

	while(true)
	{
		itrNum++;
		conv = doNextStep(itrNum,newNet);
		if (conv)
		{
			break;
		}
	}
	delete newNet;
	return 0;
}

int
Framework::init(int argc, char** argv)
{
	
	bool inputG=false;
	bool inputR=false;
	bool inputE=false;
	bool inputN=false;
	bool inputL=false;
	bool inputC=false;
	bool inputO=false;
	bool inputT=false;
	bool inputNcaRep=false;
	bool inputS=false;
	bool inputB=false;

	char gname[1025];
	char rname[1025];
	char pname[1025];
	char dname[1025];

	int optret='-';
	opterr=1;
	int oldoptind=optind;
	int condCnt=1;
	while(optret=getopt(argc,argv,"r:g:p:b:s:n:o:d:l:t:c:q:ha")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return -1;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind == 2)        //for -v 1
		{
			my_optarg=argv[oldoptind+1];	
		}
		else                             //for -v1
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'h': //Help
			{
				printHelp(argv[0]);
				return -1;
			}
			case 'o': //Output directory
			{
				strcpy(outdir,my_optarg);
				for (int ii=strlen(outdir)-1;ii>=0;ii--)
				{
					if (outdir[ii]=='/')
					{
						outdir[ii]=0;
					}
					else
					{
						break;
					}
				}
				inputO=true;
				break;
			}
			case 't': //Output directory
			{
				numThread = atoi(my_optarg);
				inputT=true;
				break;
			}
			case 'n':
			{
				numNcaRep = atoi(my_optarg);
				inputNcaRep=true;
				break;
			}
			case 's':
			{
				sparsity = atof(my_optarg);
				inputS=true;
				break;
			}
			case 'b':
			{
				beta = atof(my_optarg);
				inputB=true;
				break;
			}
			case 'l':
			{
				lambda = atof(my_optarg);
				inputL=true;
				break;
			}
			case 'c':
			{
				CV = atoi(my_optarg);
				inputC=true;
				break;
			}
			case 'q':
			{
				resIter = atoi(my_optarg);
				lShouldResume=true;
				break;
			}
			case 'a':
			{
				lResumeMERLIN=true;
				break;
			}
			case 'd':
			{
				strcpy(dname,my_optarg);
				inputE=true;
				break;
			}
			case 'p':
			{
				strcpy(pname,my_optarg);
				inputN=true;
				break;
			}
			case 'g':
			{
				strcpy(gname,my_optarg);
				inputG=true;
				break;
			}
			case 'r':
			{
				strcpy(rname,my_optarg);
				inputR=true;
				break;
			}
			default:
			{
				cerr <<"Unhandled option " << c  << endl;
				printHelp(argv[0]);
				return -1;
			}
		}
		oldoptind=optind;
	}
	if (!inputE)
	{
		cerr << "Expression matrix was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputN)
	{
		cerr << "Prior network was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputG)
	{
		cerr << "List of target genes was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputR)
	{
		cerr << "List of regulators was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputO)
	{
		cerr << "Output directory was not provided!" << endl;
		printHelp(argv[0]);
		return -1;
	}
	if (!inputB)
	{
		cerr << "Beta was not provided" << endl;
		cerr << "Setting beta to 5" << endl;
		beta = 5;
	}
	if (!inputS)
	{
		cerr << "Sparsity was not provided" << endl;
		cerr << "Setting sparsity to -5" << endl;
		sparsity = -5;
	}
	if (!inputNcaRep)
	{
		cerr << "Number of NCA replicates was not provided" << endl;
		cerr << "Setting Number of NCA replicates to 10" << endl;
		numNcaRep = 10;
	}
	if (!inputL && !inputC)
	{
		cerr << "Lambda was not provided" << endl;
		cerr << "Setting lambda to 0" << endl;
		lambda = 0;
	}
	if (inputL && inputC)
	{
		cerr << "Both Lambda and CV were provided" << endl;
		cerr << "Ignoring CV" << endl;
		CV=0;
	}
	if (!inputT)
	{
		cerr << "Number of threads was not provideed" << endl;
		cerr << "Setting to 1" << endl;
		numThread = 1;
	}

	regMngr = new VariableManager;
	regMngr->readVariables(rname);

	tgtMngr = new VariableManager;
	tgtMngr->readVariables(gname);

	readClusters(gname,moduleID);

	allGenesMngr = tgtMngr->copyMe();
	allGenesMngr->mergeVarSets(regMngr);

	evidMngr = new EvidenceManager;
	evidMngr->setVariableManager(tgtMngr);
	evidMngr->loadEvidenceFromFile_Continuous(dname);

	allEvidMngr = new EvidenceManager;
	allEvidMngr->setVariableManager(allGenesMngr);
	allEvidMngr->loadEvidenceFromFile_Continuous(dname);

	priorNet = new Graph;
	priorNet->setVariableManagers(regMngr,tgtMngr);
	priorNet->readNet(pname);

	pnet = priorNet->getNetMap();
	pnet->addSuffix("_nca");

	return 0;
}

int
Framework::writeClusters(char* aFName, map<string,int>& geneModuleID)
{
	ofstream outFile(aFName);
	for (auto itr=geneModuleID.begin();itr!=geneModuleID.end();itr++)
	{
		outFile << itr->first << "\t" << itr->second << endl;
	}
	outFile.close();
	return 0;
}

int
Framework::readClusters(char* aFName,map<string,int>& geneModuleID)
{
	geneModuleID.clear();
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
		geneModuleID[geneName]=moduleID;
	}
	inFile.close();
	return 0;
}

int
Framework::printHelp(char* name)
{
	cerr << "Usage for " << name << endl;
	cerr << "-h\t\tprint this message" << endl;
	cerr << "-d\t\tInput expression, tab delimited, [GeneName][\\t][val][\\t][val]..." << endl;
	cerr << "-r\t\tList of regulators." << endl;
	cerr << "-g\t\tList of target genes." << endl;
	cerr << "-p\t\tPrior network, [TF][\\t][TG][\\t][Confidence]" << endl;
	cerr << "-l\t\tLambda penalty, for LASSO, higher means more sparse (default 0)" << endl;
	cerr << "-c\t\tNumber of cross validation folds (instead of fixed Lambda)" << endl;
	cerr << "-s\t\tSparsity parameter, for network inference, lower means more sparse (default -5)" << endl;
	cerr << "-b\t\tBeta, importance of prior network in inference, higher means more important (default 5)" << endl;
	cerr << "-n\t\tNumber of replicates for NCA (default 10)" << endl;
	cerr << "-t\t\tNumber of threads, default 1" << endl;
	cerr << "-q\t\tContinue from an incomplete run (of NCA), input is the last complete iteration" << endl;
	cerr << "-a\t\tContinue from an incomplete run (of MERLIN) " << endl;
	cerr << "-o\t\tOutput directory, it will over write existing files" << endl;
	return 0;
}

/*
int
Framework::makeRandomSubs(int maxCnt)
{
	Matrix* d = allEvidMngr->getDataMat();
	int maxDim = d->getColCnt();
	for (int r=0;r<maxCnt;r++)
	{
		vector<int> allids;
		for (int i=0;i<maxDim;i++)
		{
			allids.push_back(i);
		}
		random_shuffle(allids.begin(),allids.end());
		vector<int> remids;
		for (int i=0;i<maxDim/2;i++)
		{
			remids.push_back(allids[i]);
		}
		allRMIDs.push_back(remids);
	}
	return 0;
}
*/

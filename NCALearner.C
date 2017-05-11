#include "NCA.H"
#include "NCALearner.H"


NCALearner::NCALearner(VariableManager* rptr, VariableManager* gptr, EvidenceManager* dptr, Graph* pptr, double l, int cv, int nr, int nt, const char* oname)
{
	lambda = l;
	CV = cv;
	numRep = nr;
	numThread = nt;
	strcpy(outdir, oname);

	regMngr = rptr->copyMe();
	tgtMngr = gptr->copyMe();

	Matrix* E = dptr->getDataMat();
	E->rowStandardize();
	expMngr = new EvidenceManager;
	expMngr->setVariableManager(tgtMngr);
	expMngr->setDataMat(E);

	Matrix* W = pptr->getAdj();
	priorNet = new Graph;
	priorNet->setVariableManagers(regMngr,tgtMngr);
	priorNet->setAdj(W);
	outregMngr=NULL;
	outtgtMngr=NULL;
	outtfaMngr=NULL;
	outnet=NULL;
	delete E;
	delete W;
}

NCALearner::~NCALearner()
{
	delete regMngr;
	delete tgtMngr;
	delete expMngr;
	delete priorNet;

	if (outregMngr!=NULL)
	{
		delete outregMngr;
	}
	if (outtgtMngr!=NULL)
	{
		delete outtgtMngr;
	}
	if (outtfaMngr!=NULL)
	{
		delete outtfaMngr;
	}
	if (outnet!=NULL)
	{
		delete outnet;
	}
}

void*
NCALearner::runOneNca(void* ptr)
{
	void** pptr = (void**)ptr;
	ThreadManager* tm = (ThreadManager*) pptr[0];

	NCA* nca = (NCA*)pptr[1];
	nca->estimate();
	VariableManager* rptr;
	VariableManager* gptr;
	Graph* aptr;
	EvidenceManager* tptr;
	nca->getData(rptr, gptr, aptr, tptr);
	delete nca;

	pthread_mutex_t& nqlock = tm->getLock();
	sem_t& nqsem = tm->getSem();

	pthread_mutex_lock(&nqlock);
	vector<VariableManager*>* allTFs   = (vector<VariableManager*>*) pptr[2];
	vector<VariableManager*>* allGenes = (vector<VariableManager*>*) pptr[3];
	vector<EvidenceManager*>* allP = (vector<EvidenceManager*>*) pptr[4];
	vector<Graph*>* allA = (vector<Graph*>*) pptr[5];
	
	allTFs->push_back(rptr);
	allGenes->push_back(gptr);
	allA->push_back(aptr);
	allP->push_back(tptr);

	pthread_mutex_unlock(&nqlock);
	sem_post(&nqsem);
	return NULL;
}

int
NCALearner::start()
{
	NCA* nca;
	vector<VariableManager*>* allTFs   = new vector<VariableManager* >;
	vector<VariableManager*>* allGenes = new vector<VariableManager* >;
	vector<EvidenceManager*>* allP = new vector<EvidenceManager*>;
	vector<Graph*>* allA = new vector<Graph*>;

	ThreadManager* thMgr = new ThreadManager(numThread,&runOneNca);
	for (int i=0;i<numRep;i++)
	{
		//cout << "Rand init " << i << endl;
		nca = new NCA(regMngr, tgtMngr, expMngr, priorNet, lambda, CV);

		void** pptr = new void*[6];
		pptr[0] = (void*) thMgr;
		pptr[1] = (void*) nca;
		pptr[2] = (void*) allTFs;
		pptr[3] = (void*) allGenes;
		pptr[4] = (void*) allP;
		pptr[5] = (void*) allA;
		thMgr->addInput((void*)pptr);
	}
	thMgr->run();
	delete thMgr;
	collapseReps(outregMngr, outtgtMngr, outtfaMngr, outnet, allTFs, allGenes, allP, allA);
	for (int i=0;i<numRep;i++)
	{
		Graph* aa = allA->at(i);
		delete aa;
		EvidenceManager* ap = allP->at(i);
		delete ap;
	}
	allA->clear();
	allP->clear();

	char mkdirCmd[1024];
	sprintf(mkdirCmd,"mkdir -p %s",outdir);
	system(mkdirCmd);
	char aoutname[1024];
	sprintf(aoutname,"%s/adj.txt",outdir);
	outnet->writeNet(aoutname);
	char poutname[1024];
	sprintf(poutname,"%s/tfa.txt",outdir);
	outtfaMngr->writeEvidence(poutname);
	return 0;
}

int 
NCALearner::countNames(vector<VariableManager*>* allVars, map<string,int>& counts)
{
	counts.clear();
	for (int i=0; i<allVars->size(); i++)
	{
		VariableManager* vMngr = allVars->at(i);
		VSET vset = vMngr->getVariableSet();
		for (auto itr=vset.begin(); itr!=vset.end(); itr++)
		{
			Variable* var = itr->second;
			string name = var->getName();
			if (counts.find(name) == counts.end())
			{
				counts[name] = 1;
			}
			else
			{
				counts[name] = counts[name]+1;
			}
		}
	}
	return 0;
}

int
NCALearner::getMissingVarIDs(VariableManager* vMngr, map<string,int>& counts, vector<int>& ids)
{
	VSET vset = vMngr->getVariableSet();
	for (auto itr=vset.begin(); itr!=vset.end(); itr++)
	{
		Variable* var = itr->second;
		string name = var->getName();
		int vID = var->getID();
		if (counts[name] < numRep)
		{
			ids.push_back(vID);
		}
	}
	return 0;
}

int
NCALearner::filterAllMissing(	
						vector<VariableManager*>* allTFs, 
						vector<VariableManager*>* allGenes, 
						vector<EvidenceManager*>* allP, 
						vector<Graph*>* allA)
{
	map<string,int> regCounts;
	map<string,int> tgtCounts;
	countNames(allTFs,regCounts);
	countNames(allGenes,tgtCounts);
	for (int i=0;i<numRep;i++)
	{
		VariableManager* rMngr = allTFs->at(i);
		VariableManager* gMngr = allGenes->at(i);
		EvidenceManager* pMngr = allP->at(i);
		Graph* net = allA->at(i);

		vector<int> regs;
		vector<int> tgts;
		getMissingVarIDs(rMngr, regCounts, regs);
		getMissingVarIDs(gMngr, tgtCounts, tgts);
		rMngr->removeVarsByID(regs);
		gMngr->removeVarsByID(tgts);
		pMngr->removeVarsByID(regs);
		net->removeRegs(regs);
		net->removeTgts(tgts);
	}
	return 0;
}

int
NCALearner::getStdCorr(Matrix* A, Matrix* B, Matrix*& C)
{
	int col = A->getColCnt();
	int row = A->getRowCnt();
	C = new Matrix(row,1);
	double v, v1, v2;
	for (int i=0;i<row; i++)
	{
		v =0;
		for (int j=0;j<col; j++)
		{
			v1 = A->getValue(i,j);
			v2 = B->getValue(i,j);
			v += (v1*v2);
		}
		v = v/((double)col);
		C->setValue(v, i, 0);
	}
	return 0;
}

int
NCALearner::collapseReps(VariableManager*& oregMngr, VariableManager*& otgtMngr, EvidenceManager*& otfaMngr, Graph*& onet,
						vector<VariableManager*>* allTFs, 
						vector<VariableManager*>* allGenes, 
						vector<EvidenceManager*>* allP, 
						vector<Graph*>* allA)
{
	filterAllMissing(allTFs, allGenes, allP, allA);

	oregMngr = allTFs->at(0)->copyMe();
	otgtMngr = allGenes->at(0)->copyMe();

	Matrix* P0 = allP->at(0)->getDataMat();
	P0->rowStandardize();
	Matrix* A0 = allA->at(0)->getAdj();

	Matrix* outP = P0->copyMe();
	Matrix* outA = A0->copyMe();

	int nRegs = P0->getRowCnt();
	int nGenes = A0->getRowCnt();
	int nSamples = P0->getColCnt();


	for (int i=1;i<numRep;i++)
	{
		Matrix* C;
		Matrix* Pi = allP->at(i)->getDataMat();
		Pi->rowStandardize();
		Matrix* Ai = allA->at(i)->getAdj();
		getStdCorr(P0, Pi, C);
		for (int j=0;j<nRegs;j++)
		{
			double crr = C->getValue(j,0);
			if (crr < 0)
			{
				//Multiply row j of P by -1
				//Multiply col j of A by -1
				for (int k=0;k<nGenes;k++)//Number of samples
				{
					double v;
					v = Ai->getValue(k,j);
					Ai->setValue(-1*v,k,j);
				}
				for (int k=0;k<nSamples;k++)//Number of samples
				{
					double v;
					v = Pi->getValue(j,k);
					Pi->setValue(-1*v,j,k);
				}
			}
		}
		outA->addWithMatrix(Ai);
		outP->addWithMatrix(Pi);
		delete Ai;
		delete Pi;
	}
	outA->multiplyScalar(1/double(numRep));
	outP->multiplyScalar(1/double(numRep));

	otfaMngr = new EvidenceManager;
	otfaMngr->setVariableManager(oregMngr);
	otfaMngr->setDataMat(outP);

	onet = new Graph;
	onet->setVariableManagers(oregMngr, otgtMngr);
	onet->setAdj(outA);

	delete A0;
	delete P0;
	delete outA;
	delete outP;

	return 0;
}

int 
NCALearner::getData(VariableManager*& rptr, EvidenceManager*& tptr, const char* suff)
{
	rptr = outregMngr->copyMeWithSuffix(suff);
	tptr = new EvidenceManager;
	tptr->setVariableManager(rptr);
	Matrix* mat = outtfaMngr->getDataMat();
	tptr->setDataMat(mat);
	delete mat;
	return 0;
}

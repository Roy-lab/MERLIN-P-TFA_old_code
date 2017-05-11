#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"


EvidenceManager::EvidenceManager()
{
	vMgr = NULL;
	dataMat = NULL;
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=0;
}

EvidenceManager::~EvidenceManager()
{
	if (dataMat!=NULL)
	{
		delete dataMat;
	}
	for (int i=0;i<evidenceSet.size(); i++)
	{
		EMAP* evidMap = evidenceSet[i];
		for (auto itr=evidMap->begin();itr!=evidMap->end(); itr++)
		{
			Evidence* evid = itr->second;
			delete evid;
		}
		evidMap->clear();
		delete evidMap;
	}
	evidenceSet.clear();
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int
EvidenceManager::removeVarsByID(vector<int>& ids)
{
	if (dataMat==NULL)
	{
		return 0;
	}
	dataMat->removeRows(ids);
	updateEvidSet();
	return 0;
}

int
EvidenceManager::readExpFromFile(const char* inname, map<string,vector<double>* >& expMap)
{
	ifstream inFile(inname);
	char buffer[10001];
	while(inFile.good())
	{
		inFile.getline(buffer,10000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		vector<double>* vals = new vector<double>;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else 
			{
				vals->push_back(stod(string(tok)));
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if (expMap.find(geneName)!=expMap.end())
		{
			cerr << "ERROR! (in readExpFromFile) Shouldn't happen, gene \"" << geneName << "\" was repeated." << endl;
		}
		expMap[geneName] = vals;
	}
	return 0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inname)
{
	if (vMgr == NULL)
	{
		cerr << "ERROR! (in loadEvidenceFromFile_Continuous): VariableManager is not set!" << endl;
		cerr << "Exiting..." << endl;
		exit(0);
	}
	map<string,vector<double>* > expMap;
	readExpFromFile(inname, expMap);
	map<int,Variable*> vset = vMgr->getVariableSet();
	vector<int> remIDs;
	for (auto itr=vset.begin(); itr!=vset.end(); itr++)
	{
		Variable* var = itr->second;
		if (expMap.find(var->getName()) == expMap.end())
		{
			remIDs.push_back(var->getID());
		}
	}
	vMgr->removeVarsByID(remIDs);
	vset = vMgr->getVariableSet();
	auto itr = expMap.begin();
	int nGenes = vset.size();
	int nSampels = itr->second->size();
	dataMat = new Matrix(nGenes,nSampels);
	for (itr=expMap.begin(); itr!=expMap.end(); itr++)
	{
		int vID = vMgr->getVarID(itr->first.c_str());
		vector<double>* vals = itr->second;
		if (vID == -1)
		{
			delete vals;
			continue;
		}
		for (int j=0;j<nSampels;j++)
		{
			dataMat->setValue(vals->at(j),vID,j);
		}
		delete vals;
	}
	expMap.clear();

	updateEvidSet();

	return Error::SUCCESS;
}

int
EvidenceManager::setDataMat(Matrix* mat)
{
	if (dataMat != NULL)
	{
		delete dataMat;
	}
	dataMat = mat->copyMe();
	updateEvidSet();
	return 0;
}

Matrix*
EvidenceManager::getDataMat()
{
	if (dataMat == NULL)
	{
		return NULL;
	}
	Matrix* temp = dataMat->copyMe();
	return temp;
}

int
EvidenceManager::updateDataMat()
{
	if (dataMat != NULL)
	{
		delete dataMat;
	}
	int nSampels = evidenceSet.size();
	int nGenes   = evidenceSet[0]->size();
	dataMat = new Matrix(nGenes,nSampels);
	for (int j=0;j<nSampels;j++)
	{
		EMAP* evidMap = evidenceSet[j];
		for (auto itr=evidMap->begin(); itr!=evidMap->end(); itr++)
		{
			int i=itr->first;
			Evidence* evid = itr->second;
			dataMat->setValue(evid->getEvidVal(),i,j);
		}
	}
	return 0;
}

int
EvidenceManager::updateEvidSet()
{
	if (dataMat == NULL)
	{
		return -1;
	}
	for (int i=0;i<evidenceSet.size(); i++)
	{
		EMAP* evidMap = evidenceSet[i];
		for (auto itr=evidMap->begin();itr!=evidMap->end(); itr++)
		{
			Evidence* evid = itr->second;
			delete evid;
		}
		evidMap->clear();
		delete evidMap;
	}
	evidenceSet.clear();
	int nSampels = dataMat->getColCnt();
	int nGenes   = dataMat->getRowCnt();
	for (int j=0;j<nSampels;j++)
	{
		EMAP* evidMap = new EMAP;
		for (int i=0;i<nGenes;i++)
		{
			Evidence* evid = new Evidence;
			evid->assocVariable(i);
			evid->setEvidVal(dataMat->getValue(i,j));
			(*evidMap)[i]=evid;
		}
		evidenceSet.push_back(evidMap);
	}
	return 0;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[trainIndex.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,trainIndex,trainIndex.size());	
		int j=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=NULL;
			if(trainIndex.find(i)!=trainIndex.end())
			{	
				evidMap=evidenceSet[randInds[j]];
				j++;
			}
			else
			{
				evidMap=evidenceSet[i];
			}
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
		string& geneName=(string&)vIter->second->getName();
		if((strcmp(geneName.c_str(),"FBgn0002631")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0000411")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0004915")==0) 
		|| (strcmp(geneName.c_str(), "FBgn0002573")==0)
		|| (strcmp(geneName.c_str(),"FBgn0005596")==0)  
		|| (strcmp(geneName.c_str(),"FBgn0035769")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0011655")==0)
		|| (strcmp(geneName.c_str(),"FBgn0000576")==0))
		{
			cout <<geneName<<"IDs";
			for(int i=0;i<trainIndex.size();i++)
			{
				cout <<"\t" <<randInds[i];
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=randEvidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=evidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

EMAP* 
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}


int
EvidenceManager::addEvidence(EMAP* evidSet)
{
	evidenceSet.push_back(evidSet);
	return 0;
}

int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::generateValidationSet(const char* vFName, int vSetSize,gsl_rng* r)
{
	ifstream inFile(vFName);
	if(inFile.good())
	{
		char buffer[256];
		while(inFile.good())
		{
			inFile.getline(buffer,255);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			int dId=atoi(buffer);
			validationIndex[dId]=0;
		}
		inFile.close();
	}
	else
	{
		populateRandIntegers(r,validationIndex,evidenceSet.size(),vSetSize);
		ofstream oFile(vFName);
		for(INTINTMAP_ITER vIter=validationIndex.begin();vIter!=validationIndex.end();vIter++)
		{
			oFile << vIter->first << endl;
		}
		oFile.close();
	}
	return 0;
}


int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::splitData(int s)
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	if(foldCnt==1)
	{
		testStartIndex=-1;
		testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		randseed=getpid();
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());	
		gsl_rng_free(r);
		//cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[eInd]=0;
		}
		else
		{
			trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}

INTINTMAP& 
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}

INTINTMAP& 
EvidenceManager::getTestSet()
{
	return testIndex;
}

INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}


int 
EvidenceManager::standardizeData()
{
	VSET& varSet=vMgr->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		double mean=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			mean=mean+(*evidMap)[vIter->first]->getEvidVal();
		}
		mean=mean/evidenceSet.size();
		double std=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			double diff=mean-(*evidMap)[vIter->first]->getEvidVal();
			std=std+(diff*diff);
		}
		std=sqrt(std/(evidenceSet.size()-1));
		//Standardize
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			double tval=evid->getEvidVal();
			double sval=(tval-mean)/std;
			evid->setEvidVal(sval);
		}
	}
	return 0;
}

int
EvidenceManager::partitionData(int numberOfComponents,map<int,EvidenceManager*>& evMgrSet,int& rseed,map<int,INTINTMAP*>& datasetInd)
{
	int datasubsetSize=evidenceSet.size()/numberOfComponents;
	int* randInds=new int[evidenceSet.size()];
	//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//randseed=getpid();
	randseed=18721;
	rseed=randseed;
	gsl_rng_set(r,randseed);
	populateRandIntegers(r,randInds,evidenceSet.size());	
	//for(int i=0;i<evidenceSet.size();i++)
	//{
	//	cout << randInds[i] << endl;
	//}
	gsl_rng_free(r);
	int ind=0;
	for(int n=1;n<=numberOfComponents;n++)
	{
		int startIndex=(n-1)*(datasubsetSize);
		int endIndex=n*datasubsetSize;
		if(n==numberOfComponents)
		{
			endIndex=evidenceSet.size();
		}
		EvidenceManager* localManager=new EvidenceManager;
		INTINTMAP* origIDs=new INTINTMAP;
		int currind=(int)pow(2.0,ind);
		datasetInd[currind]=origIDs;
		evMgrSet[currind]=localManager;
		int dId=0;
		for(int e=startIndex;e<endIndex;e++)
		{
			//int rId=randInds[e];
			int rId=e;
			EMAP* evidSet=evidenceSet[rId];
			localManager->addEvidence(evidSet);
			(*origIDs)[dId]=rId;
			dId++;
		}
		//localManager->updateDataMat();
		ind++;
	}
	delete[] randInds;
	return 0;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	map<int,int> temp;
	INTINTMAP_ITER tIter=populateFrom.begin();
	for(int i=0;i<size;i++)
	{
		int tid=tIter->first;
		temp[i]=tid;
		tIter++;
	}
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=temp[rind];
	}
	usedInit.clear();
	return 0;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, INTINTMAP& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(randInds.find(rind)!=randInds.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		randInds[rind]=0;
	}
	return 0;
}

int 
EvidenceManager::populateRandIntegers(gsl_rng* r, vector<int>& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds.push_back(rind);
	}
	return 0;
}

int
EvidenceManager::writeEvidence(const char* oname)
{
	ofstream oFile(oname);
	int row = dataMat->getRowCnt();
	int col = dataMat->getColCnt();
	for (int i=0;i<row;i++)
	{
		Variable* var = vMgr->getVariableAt(i);
		oFile << var->getName();
		for (int j=0;j<col;j++)
		{
			oFile << "\t" << dataMat->getValue(i,j);
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

//It assumes that VariableManager is already updated
int 
EvidenceManager::updateEvidence(EvidenceManager* eMngr)
{
	int nGenes = dataMat->getRowCnt();
	int nSampels = dataMat->getColCnt();
	//Resize the dataMat if the size of VariableManager and dataMat don't match
	if (nGenes < vMgr->getVarCnt())
	{
		int oldNGenes = nGenes;
		nGenes = vMgr->getVarCnt();
		Matrix* temp = new Matrix(nGenes,nSampels);
		for (int i=0;i<oldNGenes;i++)
		{
			for (int j=0;j<nSampels;j++)
			{
				temp->setValue(dataMat->getValue(i,j),i,j);
			}
		}
		delete dataMat;
		dataMat = temp;
	}
	VSET their_vset = eMngr->vMgr->getVariableSet();
	for (auto itr=their_vset.begin(); itr!=their_vset.end(); itr++)
	{
		Variable* their_var = itr->second;
		int their_i = their_var->getID();
		string name = their_var->getName();
		int our_i   = vMgr->getVarID(name);
		if (our_i == -1)
		{
			cerr << "ERROR! (updateEvidence), two VariableManagers don't match!" << endl;
			cerr << "Exiting..." << endl;
			exit(0);
		}
		for (int j=0;j<nSampels;j++)
		{
			dataMat->setValue(eMngr->dataMat->getValue(their_i,j),our_i,j);
		}
	}
	updateEvidSet();
	return 0;
}

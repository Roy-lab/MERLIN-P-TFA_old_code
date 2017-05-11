#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include <gsl/gsl_randist.h>

#include "Error.H"
#include "Evidence.H"
#include "VariableManager.H"
#include "Matrix.H"
#include "CommonTypes.H"

using namespace std;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
typedef map<int,Evidence*> EMAP;
typedef map<int,Evidence*>::iterator EMAP_ITER;

typedef vector<EMAP*> EVIDENCE_SET;

class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();

		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);
		int randomizeEvidence(gsl_rng*);
		int getNumberOfEvidences();
		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);
		EMAP* getRandomEvidenceAt(int);

		int setFoldCnt(int);
		int generateValidationSet(const char*,int vSetSize,gsl_rng* r);
		int setPreRandomizeSplit();
		int splitData(int);
		INTINTMAP& getTrainingSet();
		INTINTMAP& getTestSet();
		INTINTMAP& getValidationSet();
		int populateRandIntegers(gsl_rng*, INTINTMAP&,int,int);
		int populateRandIntegers(gsl_rng*, vector<int>&,int,int);
		int standardizeData();
		int partitionData(int,map<int,EvidenceManager*>&,int&,map<int,INTINTMAP*>&);

		int setDataMat(Matrix* mat);
		Matrix* getDataMat();

		int removeVarsByID(vector<int>& ids);

		int writeEvidence(const char* oname);
		int updateEvidence(EvidenceManager* eMngr);
	private:
		int readExpFromFile(const char* inname, map<string,vector<double>* >& expMap);
		//When someone updates Matrix, we update the evidenceSet
		int updateEvidSet();
		//When someone updates the evidenceSet, we update Matrix
		int updateDataMat();
		//Only EvidenceManager can call this, so we can make sure to call update afterward
		int addEvidence(EMAP*);

		int populateRandIntegers(gsl_rng*, int*,int);
		int populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size);
		EVIDENCE_SET evidenceSet;
		EVIDENCE_SET randEvidenceSet;
		VariableManager* vMgr;
		INTINTMAP trainIndex;
		INTINTMAP testIndex;
		INTINTMAP validationIndex;
		int foldCnt;
		bool preRandomizeSplit;
		int randseed;

		Matrix* dataMat;
};

#endif //_EVIDENCE_MANAGER

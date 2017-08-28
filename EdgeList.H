#ifndef _EDGELIST
#define _EDGELIST

#include "CommonTypes.H"
#include <vector>
#include <map>
#include <string>
using namespace std;

class EdgeList
{
	public:
		EdgeList();
		~EdgeList();
		int addEdge(string& tname, string& gname, double v);
		NET_T* getNet();
		int setNet(NET_T*);
		int countEdges();
		int setConsensus(vector<EdgeList*>*);
		EdgeList* compyMe();
		int getOverlap(EdgeList* onet);
		int addSuffix(const char* suff);
		int setUnion(vector<EdgeList*>* outnets);
		int writeNet(char* oname);
		int readNet(const char* inname);
		EdgeList* removeSuffix(const char* suff);
		EdgeList* copyMe();
		EdgeList* percentileRank();
	private:
		int deleteNet();
		map<string,map<string,double>* >* net;
		struct customLess_T
		{
			bool operator()(pair<string,double> a, pair<string,double> b) const
			{   
				return a.second > b.second;
			}   
		};
		int rankeOneTF(map<string,double>* elist, map<string,double>* nelist);
};

#endif //_EDGELIST

#include <vector>

bool any(bool *G,int p);
void initialiseGraph(bool* G,int p);
int getNextRowWithConnections(int startrow,bool* G,int p);
void getRowConnections(int row,bool* G,int p,int* connections);
int getNextRowWithConnections(int startrow,bool* G,int p);
vector<int> getOtherConnections(int j,int* connections,int p);
vector<int> getSeqVector(int ord);
vector<int> getSubset(std::vector<int> set,std::vector<int> subsetind);


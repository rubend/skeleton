bool any(bool *G,int p);
void initialiseGraph(bool* G,int p);
int getNextRowWithConnections(int startrow,bool* G,int p);
void getRowConnections(int row,bool* G,int p,int* connections);
int getNextRowWithConnections(int startrow,bool* G,int p);
std::vector<int> getOtherConnections(int j,int* connections,int p);
std::vector<int> getSeqVector(int ord);
std::vector<int> getSubset(std::vector<int> set,std::vector<int> subsetind);


int it = as<int> (i) -1;
int jt = as<int> (j) -1;

//R vector convert to std::vector k
IntegerVector l(k);
std::vector<int> kt(l.size());

for (int i = 0; i < l.size(); ++i)
  {
    kt[i] = l[i] -1;
    
  }
NumericMatrix Corr(C);

cout << "own pcor value " << pcorOrder(it,jt,kt,Corr) << endl;


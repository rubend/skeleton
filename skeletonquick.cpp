// now this main will be implemented as the skeleton function! :p
  
int p = as<int>(pt);
double alpha = as<double>(alphat);
int m_max = as<int>(m_maxt);
NumericMatrix Corr(C);
long n = as<long>(nt);
  
//using c++ datatypes and trying to write all functions myself, 
//i.e. not calling R from c++
//how to store graph? Use matrix?
bool G[p*p];// only store upper triangle? -> store as vector?
  
initialiseGraph(G,p);

int row;// the current row we are studying the connections of
int connections[p-1]; // stores the connections the current row has (maximum p-1)
//static allocation ==> efficient
int sizeothers,x,y;
double pval;
std::vector<int> subset(0,0);//declare in beginning and increase size every loop run
std::vector<int> others(0);
//ord = 0 gives problems with subset size and the like, right? is it really needed?
int iter =0;

//SPEEDUP: parallelise with openmp!
for (int ord = 0; ord <= m_max; ++ord)
  //additional check for if any connections are left might be useful, especially if mmax is huge
  {
    cout << "ord = " << ord << endl;
    
    //look for next row with connections then iterate over the remaining connections
    //alternative: save all connections explicitly in double array and loop over those.
    row = 0; // reset to zero
    
    row = getNextRowWithConnections(row,G,p); // row == x
   //cout << "row = " << row << endl;

    //exit(0);
    
    
      
    while(row != -1)
      {
	x = row; // just temporary variable for clarity when calling pcororder
	getRowConnections(row,G,p,connections); // getting the connections belonging to the row
	for (int i = 0; i < p-1; ++i)
	  {
	    //one of the remaining edge tests
	    y=connections[i];
	      
	    //so now we are gonna check the correlation between row and connections[i]
	    //in respect to every other subset of lengths ord of the remaining connections
	    if(y == -1) break;//reached end of connections
	    //y == connections[i]
	    getOtherConnections(&others,i,connections,p);
	    sizeothers = others.size();

	    if (sizeothers < ord)
	      {
		continue;//goto next loop iteration
	      }
	      
	    //initial subset, TODO is there a more efficient way? Builtin way for this?
	      
	    getSeqVector(&subset,ord);//get sequence vector of size ord in subset
	    
	    //problem is that subset has size 0, getting the 0th element 
	    //--> gives shit! --> crappy way of doing stuff.
	    
	    //temporary fix, might be a cleaner way to do this
	    
	    while(subset.size() == 0 || subset[0] != -1)
	      {
		//subset.size() == 0 is possible if ord == 0, subset[0] == -1 is the stop condition
		
		std::vector<int> k = getSubset(others,subset); //does this work? check! otherwise write own function that does this.
		//pval = pcorOrder(x,y,k,C)
		//cout << "k " << k.back() << endl;
		
		//WRONG! pval = pcorOrder(x,y,k,Corr);
		pval = gaussCItest(x,y,k,Corr,n);
		
		if(k.size()>0)
		  {
		    cout <<"k ="<< k[0]<<endl;
		    cout << "x = "<<x+1 <<"y = " << y+1 << " S[0]=" << k[0] << " : pval = " << pval << endl;
		  }else
		  {
		    
		    cout << "x = "<<x+1 <<"y = " << y+1 << " S[0]=" << " : pval = " << pval << endl;
		  }
		
		cout << "pval = " << pval << endl;
		
		iter++;
		if (x == 2 && y == 3 && k.size()==1 )
		  {
		    cout << "pval in case x 2 y 3 k " << k[0] <<" =" << pval << endl;
		  }
		//if ((ord == 0 && pval <  alpha) || (ord >0 && pval >= alpha))
		if (pval >=  alpha)
		  {
		    cout << "x " << x << " and y " << y << endl;
		    
		    cout << "FALSE" << endl;
		    
		    //independent
		    G[p*x+y]=false;G[p*y+x]=false;
	
		    break; //no more checking to be done
		  }
		subset = getNextSet(sizeothers,ord,subset);  
	      }
	  }
	
	row = getNextRowWithConnections(row+1,G,p);
      }
  }

return wrap(convertToLogical(G,p)); // return graph matrix, in later stage return a more complete object
